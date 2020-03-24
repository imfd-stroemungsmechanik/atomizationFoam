/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dropletCloud.H"
#include "collisionModel.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::collisionModel::collideDroplets
(
    const scalar dt,
    droplet& p1,
    droplet& p2,
    scalar& m1,
    scalar& m2
)
{
    bool coalescence = false;

    const vector pos1(p1.position());
    const vector pos2(p2.position());

    const vector& U1 = p1.U();
    const vector& U2 = p2.U();

    vector URel(U1 - U2);

    vector d(pos2 - pos1);
    scalar magd = mag(d);

    scalar vAlign = URel & (d/(magd + ROOTVSMALL));

    // Droplets travel towards each other
    if (vAlign > 0)
    {
        const scalar d1 = p1.d();
        const scalar d2 = p2.d();

        scalar sumD = d1 + d2;

        if (vAlign*dt > magd - 0.5*sumD)
        {
            scalar magU1 = mag(U1) + ROOTVSMALL;
            scalar magU2 = mag(U2) + ROOTVSMALL;
            vector n1 = U1/magU1;
            vector n2 = U2/magU2;

            scalar n1n2 = n1 & n2;
            scalar n1d = n1 & d;
            scalar n2d = n2 & d;

            scalar det = 1.0 - sqr(n1n2);

            scalar alpha = GREAT;
            scalar beta = GREAT;

            if (mag(det) > 1.0e-4)
            {
                beta = -(n2d - n1n2*n1d)/det;
                alpha = n1d + n1n2*beta;
            }

            alpha /= magU1*dt;
            beta /= magU2*dt;

            // Check if collision is possible within this timestep
            if ((alpha > 0) && (alpha < 1.0) && (beta > 0) && (beta < 1.0))
            {
                vector p1c = pos1 + alpha*U1*dt;
                vector p2c = pos2 + beta*U2*dt;

                scalar closestDist = mag(p1c - p2c);

                scalar collProb =
                    pow(0.5*sumD/max(0.5*sumD, closestDist), cSpace_)
                   *exp(-cTime_*mag(alpha - beta));
                scalar prob = ranGen_.sample01<scalar>();

		// collision occurs
                if (prob < collProb)
                {
                    if (d1 > d2)
                    {
                        coalescence = this->collideSorted(dt, p1, p2, m1, m2);
                    }
                    else
                    {
                        coalescence = this->collideSorted(dt, p2, p1, m2, m1);
                    }
                }
            }
        }
    }

    return coalescence;
}


bool Foam::collisionModel::collideSorted
(
    const scalar dt,
    droplet& p1,
    droplet& p2,
    scalar& m1,
    scalar& m2
)
{
    const scalar d1 = p1.d();
    const scalar d2 = p2.d();

    const vector& U1 = p1.U();
    const vector& U2 = p2.U();

    vector URel = U1 - U2;
    scalar magURel = mag(URel);

    scalar mTot = m1 + m2;

    scalar gamma = d1/max(ROOTVSMALL, d2);
    scalar f = pow3(gamma) + 2.7*gamma - 2.4*sqr(gamma);

    scalar dAve = sqrt(d1*d2);
    scalar WeColl = 0.5*rhop_*sqr(magURel)*dAve/max(ROOTVSMALL, sigma_);

    scalar coalesceProb = min(1.0, 2.4*f/max(ROOTVSMALL, WeColl));
    scalar prob = ranGen_.sample01<scalar>();

    // Coalescence
    if (prob < coalesceProb)
    {
        p1.U() = (m1*U1 + m2*U2)/(m1+m2);
        
        m1 += m2;
        m2 = -1;

        return true;
    }
    // Stretching separation
    else
    {
        scalar gf = sqrt(prob) - sqrt(coalesceProb);
        scalar denom = 1.0 - sqrt(coalesceProb);
        if (denom < 1.0e-5)
        {
            denom = 1.0;
        }
        gf /= denom;

        // If gf negative, this means that coalescence is turned off
        // and these parcels should have coalesced
        gf = max(0.0, gf);

        vector mr = m1*U1 + m2*U2;
        vector v1p = (mr + m2*gf*URel)/mTot;
        vector v2p = (mr - m1*gf*URel)/mTot;

        p1.U() = v1p;
        p2.U() = v2p;        

	return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::collisionModel::collisionModel
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_
    (
        IOobject
        (
            "cloudProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    active_(readBool(dict_.lookup("collision"))),
    ranGen_(clock::getTime() + pid()),
    cTime_(dict_.lookupOrDefault("cTime", 1.0)),
    cSpace_(dict_.lookupOrDefault("cSpace", 0.3)),
    rhop_(0.0),
    sigma_(0.0)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::collisionModel::update
(
    dropletCloud& cloud,
    const scalar dt
)
{
    if (!active_)
    {
        return;
    }

    // Store cloud droplet properties
    rhop_ = cloud.rhop();
    sigma_ = cloud.sigma();

    // Create pointer to cloud
    dropletCloud *ptrCloud = const_cast<dropletCloud*>(&cloud);    

    // Create the occupancy list for the cells
    labelList occupancy(mesh_.nCells(), 0);
    forAllIter(Cloud<droplet>, *ptrCloud,  iter)
    {
        occupancy[iter().cell()]++;
    }

    // Initialize the sizes of the lists of parcels in each cell
    CompactListList<droplet*> pInCell(occupancy);

    // Reset the occupancy to use as a counter
    occupancy = 0;

    // Set the parcel pointer lists for each cell
    forAllIter(Cloud<droplet>, *ptrCloud, iter)
    {
        pInCell(iter().cell(), occupancy[iter().cell()]++) = &iter();
    }

    // Perform droplet-droplet collision per cell
    for (label celli=0; celli < mesh_.nCells(); celli++)
    {
        UList<droplet*> pInCelli(pInCell[celli]);

        if (pInCelli.size() >= 2)
        {
            forAll(pInCelli, i)
            {
                for (label j=i+1; j<pInCelli.size(); j++)
                {
                    droplet& p1 = *pInCelli[i];
                    droplet& p2 = *pInCelli[j];

                    scalar m1 = rhop_*(4.0/3.0)*constant::mathematical::pi*pow(p1.d()/2.0, 3.0);
                    scalar m2 = rhop_*(4.0/3.0)*constant::mathematical::pi*pow(p2.d()/2.0, 3.0);

                    bool massChanged = collideDroplets(dt, p1, p2, m1, m2);
    
                    if (massChanged)
                    {
                        if (m1 > ROOTVSMALL)
                        {
                            p1.d() = cbrt(6.0*m1/(rhop_*constant::mathematical::pi));
                        }
                        else
                        {
                            p1.d() = -1;
                        }
                        if (m2 > ROOTVSMALL)
                        {
                            p2.d() = cbrt(6.0*m2/(rhop_*constant::mathematical::pi));
                        }
                        else
                        {
                            p2.d() = -1;
                        }
                    }
                }
            }
        }
    }

    // Delete all droplets with negative diameter
    forAllIter(Cloud<droplet>, *ptrCloud, iter)
    {
        droplet& p = iter();
        
        if (p.d() < VSMALL)
        {
            cloud.deleteParticle(p);
        }
    }
}

// ************************************************************************* //
