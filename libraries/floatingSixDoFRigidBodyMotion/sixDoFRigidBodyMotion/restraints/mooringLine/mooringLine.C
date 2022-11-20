/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "mooringLine.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(mooringLine, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        mooringLine,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::mooringLine::mooringLine
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict),
    anchor_(),
    refAttachmentPt_(),
    mass_(),
    length_(),
    gMag_(),
    gravityVector_(),
    unitVert_(),

    nCells_(101),

    simpleState_("simpleState"),
    restingState_("restingState"),
    hangingState_("hangingState"),
    mooringState_("NULL")
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::mooringLine::~mooringLine()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::pointField Foam::sixDoFRigidBodyMotionRestraints::mooringLine::centreLine
(
    const point& start,
    const point& end
) const
{
    setState(start, end);

    Info << "Writing for state: " << mooringState_ << endl;

    pointField cl(nCells_ + 1, vector::zero);

    if (mooringState_ == simpleState_)
    {
        // Create the catenary shape
        catenaryShape cs(start, end, length_, mass_, gravityVector_);

        // Calculate the centreline points
        cs.centreLine(cl);
    }
    else if (mooringState_ == restingState_)
    {
        // Find the touchdown point
        point E = end;
        scalar newLength(0);
        restingLength(start, E, newLength);

        // Calculate the catenary shape between fairlead and seabed
        catenaryShape cs(start, E, newLength, mass_, gravityVector_);

        //  Remove one point from the centreline and calculate centreline points
        cl.setSize(nCells_);

        cs.centreLine(cl);

        // Add the endpoint to the centreline - the horizontal part of the line
        cl.setSize(nCells_ + 1);
        cl[nCells_] = end;
    }
    else if (mooringState_ == hangingState_)
    {
        point midPoint = start + ((end - start) & unitVert_)*unitVert_;

        label N0 = nCells_/2;
        label N1 = nCells_ + 1 - N0;

        for (int i = 0; i < N0; i++)
        {
            scalar factor = static_cast<scalar>(i)/static_cast<scalar>(N0);
            cl[i] = start + factor*(midPoint - start);
        }

        for (int i = 0; i < N1; i++)
        {
            scalar factor = static_cast<scalar>(i)/static_cast<scalar>(N1 - 1);
            cl[i + N0] = midPoint + factor*(end - midPoint);
        }
    }
    else
    {
        notImplemented("The state is not defined.");
    }


    // Expand the centerline to a four points with the centre line in the middle
    pointField pp(cl.size()*4, point::zero);

    // Compute the points on the catenary mesh
    vector normal1(vector::zero);
    vector normal2(vector::zero);

    label pCount = 0;

    label N = cl.size();

    // Thickness of the mesh
    scalar thickness_(readScalar(sDoFRBMRCoeffs_.lookup("thickness")));

    for (int j = 0; j < N - 1; j++)
    {
        vector unit(unitVert_);

         vector vec = cl[j + 1] - cl[j];
           scalar vecMag = Foam::mag(vec);
          vec /= vecMag;

          if (Foam::mag(Foam::mag(vec & unit) - Foam::mag(vec)) < SMALL)
          {
              unit += vector(1, 0, 0);
              unit /= Foam::mag(unit);
          }

           normal1 = unit - (vec & unit)*vec;
        normal1 /= Foam::mag(normal1);
        normal2 = (normal1 ^ vec);
        normal2 /= Foam::mag(normal2);

        normal1 *= thickness_*0.5;
        normal2 *= thickness_*0.5;

        // ++
        pp[pCount++] = cl[j] + 0.5*normal1 + 0.5*normal2;

        // +-
        pp[pCount++] = cl[j] + 0.5*normal1 - 0.5*normal2;

        // --
        pp[pCount++] = cl[j] - 0.5*normal1 - 0.5*normal2;

        // -+
        pp[pCount++] = cl[j] - 0.5*normal1 + 0.5*normal2;
    }

    // Adding the last point and reuse the normals
    // ++
    pp[pCount++] = cl[N - 1] + 0.5*normal1 + 0.5*normal2;

    // +-
    pp[pCount++] = cl[N - 1] + 0.5*normal1 - 0.5*normal2;

    // --
    pp[pCount++] = cl[N - 1] - 0.5*normal1 - 0.5*normal2;

    // -+
    pp[pCount++] = cl[N - 1] - 0.5*normal1 + 0.5*normal2;

    return pp;
}


void Foam::sixDoFRigidBodyMotionRestraints::mooringLine::setState
(
    const point& start,
    const point& end
) const
{
    vector span = (start - end);
    span -= (span & unitVert_)*unitVert_;

    // Making sure that fairlead is not above anchor point
    if (Foam::mag(span) < 0.001*length_)
    {
        mooringState_ = hangingState_;

        return;
    }

    // Investigate, whether the angle at the fairlead exceeds 88 degrees.
    // Inside scope for cleanness
    {
        // Maximum top slope is set to 88 degrees
        scalar topSlope(Foam::tan(maxAngle_/180.*M_PI));
        scalar lk2 = Foam::asinh(topSlope);

        // Shape parameter (assuming 88 degrees at the top)
        scalar k =
            1.0/Foam::mag((start - end) & unitVert_)*(Foam::cosh(lk2) - 1.0);

        // Half length of the sag
        scalar L2 = 1/k*topSlope;

        // Half span distance of the sag
        scalar l2 = Foam::asinh(L2*k)/k;

        scalar maximumLength = L2 + Foam::mag(span) - l2;

        if (maximumLength < length_)
        {
            mooringState_ = hangingState_;

            return;
        }
    }

    // Investigate, whether the "top" point is in between the fairlead and the
    // anchor points.
    // It is now 'safe' to create the catenary shape, because it is known that
    // the line is not in a hanging state
    catenaryShape cs(start, end, length_, mass_, gravityVector_);

    if (cs.isUShaped())
    {
        mooringState_ = restingState_;
    }
    else
    {
        mooringState_ = simpleState_;
    }
}


void Foam::sixDoFRigidBodyMotionRestraints::mooringLine::restingLength
(
    const point& start,
    point& end,
    scalar& lNew
) const
{
    scalar lengthNew = 0;

    // Horizontal and vertical span
    vector span = end - start;
    span -= (span & unitVert_)*unitVert_;

    scalar h = Foam::mag((end - start) & unitVert_);

    // Define limits for the search
    scalar alphaMax = maxAngle_/180.0*M_PI;

    // The lower limit is based on the angle for the catenary shape, that has
    // a horizontal gradient at the anchor point
    scalar alphaMin = 0.0;

    // Enclosed in a scope for cleanness
    {
        scalar kMin = 1.0e-10;
        // The factor 1.1 is just to be sure
        scalar kMax = 1.1*Foam::mag(Foam::asinh(Foam::tan(alphaMax)))/Foam::mag(span);

        while (true)
        {
            scalar kAv = 0.5*(kMin + kMax);

            scalar valAv = -h*kAv - (1.0 - Foam::cosh(kAv*Foam::mag(span)));

            if (Foam::mag(valAv) < 1e-7)
            {
                alphaMin = Foam::atan(Foam::sinh(kAv*Foam::mag(span)));
                break;
            }

            if (valAv < 0.0)
            {
                kMin = kAv;
            }
            else
            {
                kMax = kAv;
            }
        }
    }

    // Find the resting length and the touchdown point
    label count = 0;

    while (true)
    {
        scalar alphaAv = 0.5*(alphaMin + alphaMax);

        scalar slope = tan(alphaAv);
        scalar lk2 = Foam::asinh(slope);
        scalar k = 1.0/h*(Foam::cosh(lk2) - 1);
        scalar L2 = 1.0/(k)*slope;
        scalar l2 = asinh(L2*k)/k;

        scalar totalLength = L2 + Foam::mag(span) - l2;

        if (totalLength < length_)
        {
            alphaMin = alphaAv;
        }
        else
        {
            alphaMax = alphaAv;
        }

        if (Foam::mag(totalLength - length_) < 1.0e-7*length_ || count == 100)
        {
            lengthNew = L2;

            end = start + l2*span/Foam::mag(span) - h*unitVert_;

            break;
        }
        count++;
    }

    lNew = lengthNew;
}


void Foam::sixDoFRigidBodyMotionRestraints::mooringLine::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    // Set restraint position, force and moment
    restraintPosition = motion.transform(refAttachmentPt_);
    restraintForce = vector::zero;
    restraintMoment = vector::zero;

    // Rename
    point pos0 = restraintPosition;
    point pos1 = anchor_;

    // Set the state of the mooring line
    setState(restraintPosition, anchor_);

    // Compute the force as a function of the state
    if (mooringState_ == simpleState_)
    {
        // Get access to the catenary shape
        catenaryShape cs(pos0, pos1, length_, mass_, gravityVector_);

        // Get the forces
        vector H0 = cs.H0();
        vector R0 = cs.R0();

        restraintForce = (H0 + R0);
    }
    else if (mooringState_ == restingState_)
    {
        // Find the touchdown point
        point E = pos1;
        scalar newLength(0);
        restingLength(pos0, E, newLength);

        // Calculate the catenary shape between fairlead and seabed
        catenaryShape cs(pos0, E, newLength, mass_, gravityVector_);

        // Get the forces
        vector H0 = cs.H0();
        vector R0 = cs.R0();

        restraintForce = (H0 + R0);
    }
    else if (mooringState_ == hangingState_)
    {
        // The vertical force is equal to the submerged weight of the line
        // hanging from the fairlead to the seabed (assumed at the level of
        // pos1)
        restraintForce = ((pos1 - pos0) & unitVert_)*unitVert_*mass_*gMag_;
    }
    else
    {
        notImplemented("The force calculation is not implemented");
    }

    if (motion.report())
    {
        Info<< " state: " << mooringState_
            << " force " << restraintForce
            << " position " << restraintPosition
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::mooringLine::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.lookup("anchor") >> anchor_;
    sDoFRBMRCoeffs_.lookup("refAttachmentPt") >> refAttachmentPt_;
    sDoFRBMRCoeffs_.lookup("massPerLength") >> mass_;
    sDoFRBMRCoeffs_.lookup("lineLength") >> length_;
    sDoFRBMRCoeffs_.lookup("gravityVector") >> gravityVector_;

    gMag_ = Foam::mag(gravityVector_);
    unitVert_ = gravityVector_/gMag_;
    unitVert_ = Foam::cmptMultiply(unitVert_, unitVert_);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::mooringLine::write
(
    Ostream& os
) const
{
    os.writeKeyword("anchor")
        << anchor_ << token::END_STATEMENT << nl;

    os.writeKeyword("refAttachmentPt")
        << refAttachmentPt_ << token::END_STATEMENT << nl;

    os.writeKeyword("massPerLength")
        << mass_ << token::END_STATEMENT << nl;

    os.writeKeyword("lineLength")
        << length_ << token::END_STATEMENT << nl;

    os.writeKeyword("gravityVector")
        << gravityVector_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
