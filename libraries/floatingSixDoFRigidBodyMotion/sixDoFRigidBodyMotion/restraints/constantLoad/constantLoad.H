/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::sixDoFRigidBodyMotionRestraints::constantLoad

Description
    sixDoFRigidBodyMotionRestraints model.  Constant load (force and torque).

    Example specification in dynamicMeshDict:
    \verbatim
     restraints
     {
        exampleLoad
        {
            sixDoFRigidBodyMotionRestraint constantLoad;
            applicationPt   (0 0 0); // Force application point
            movePt          true; // Application point moves with rigid body
            force           (1 2 3); // Constant force (global frame)
            torque          (-1 -2 -3); // Constant moment (global frame)  
        }

        weightLoad
        {
            sixDoFRigidBodyMotionRestraint constantLoad;
            applicationPt   (0 0 0); // Initial CoG
            movePt          true; // Move CoG with rigid body
            force           (0 -100 0); // Weight (global frame)
            torque          (0 0 0); // None  
        }
     }
    \endverbatim

SourceFiles
    constantLoad.C

Author
    Pere Frontera Pericàs, TUDelft
\*---------------------------------------------------------------------------*/

#ifndef constantLoad_H
#define constantLoad_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class constantLoad Declaration
\*---------------------------------------------------------------------------*/

class constantLoad
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data

        //- Point where the constant force is being applied
        point applicationPt_;

        //- Force vector
        vector force_;

        //- Torque vector
        vector torque_;

        //- Switch to move the application point with the rigid body
        bool movePt_;

public:

    //- Runtime type information
    TypeName("constantLoad");


    // Constructors

        //- Construct from components
        constantLoad
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new constantLoad(*this)
            );
        }


    //- Destructor
    virtual ~constantLoad();


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
