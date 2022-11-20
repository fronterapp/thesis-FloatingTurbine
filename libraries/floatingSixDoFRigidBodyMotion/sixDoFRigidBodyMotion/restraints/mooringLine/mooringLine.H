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

Class
    Foam::sixDoFRigidBodyMotionRestraints::mooringLine

Description
    sixDoFRigidBodyMotionRestraints model.  Mooring Line

    Example specification in dynamicMeshDict:
    \verbatim
     restraints
     {
        catenaryLine
        {
            sixDoFRigidBodyMotionRestraint mooringLine;
            anchor              (0 1 -1); // Anchor point, where the line is attached to an immovable object
            refAttachmentPt     (0 0 0); // Reference point of attachment to the solid bodyy
            massPerLength       1; // Mass per unit length
            lineLength          10; // Length of the mooring line
            gravityVector       (0 0 -9.81); // Gravity vector
            thickness           1; // Line thickness, used only for visualization purposes
        }
     }
    \endverbatim

SourceFiles
    mooringLine.C

Author
    Niels G. Jacobsen, Deltares

\*---------------------------------------------------------------------------*/

#ifndef mooringLine_H
#define mooringLine_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "fvMesh.H"
#include "catenaryShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class mooringLine Declaration
\*---------------------------------------------------------------------------*/

class mooringLine
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data

        //- Anchor point, where the spring is attached to an immovable
        //  object
        point anchor_;

        //- Reference point of attachment to the solid body
        point refAttachmentPt_;

        //- Mass per unit length
        scalar mass_;

        //- Length of the mooring line
        scalar length_;

        //- Magnitude of the acceleration due to gravity
        scalar gMag_;

        //- Gravity vector
        vector gravityVector_;

        //- Unit vertical direction
        vector unitVert_;

        label nCells_;

        //- Maximum slope of the line with respect to vertical
        const scalar maxAngle_ = 88.0;

        //- States
        word simpleState_;
        word restingState_;
        word hangingState_;

        mutable word mooringState_;

        //- Checks the state of the line
        void setState(const point&, const point&) const;

        //- Returns the length of suspended line in restingState_
        void restingLength(const point&, point&, scalar&) const;

public:

    //- Runtime type information
    TypeName("mooringLine");


    // Constructors

        //- Construct from components
        mooringLine
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new mooringLine(*this)
            );
        }


    //- Destructor
    virtual ~mooringLine();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;

        //- Write the centre line
        virtual pointField centreLine(const point&, const point&) const;

        void meshDimensions
        (
            label& nPoints,
            label& nFaces,
            label& nBndFaces,
            label& nCells
        ) const
        {
            // Defined from the individual restraints
            nCells = nCells_;

            // Always the same
            nPoints   = 4*(nCells_ + 1);
            nFaces    = 4*nCells_ + nCells_ + 1;
            nBndFaces = 4*nCells_ + 2;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
