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
    Foam::sixDoFRigidBodyMotionRestraints::gyroscopicMoment


Description
    sixDoFRigidBodyMotionRestraints model. 
    Accounts for the gyroscopic effect imparted by a rotating body which 
    is attatched to the six-DoF rigid body, whose rotation is not modelled.

    For floating turbines, this restraint includes the gyroscopic moment caused
    by the rotating turbine since the latter is not included in the rigid body
    definition.

    The present approach is based on the assumtption that the rotation motion
    takes place at constant angular speed. 

    Reference:
    \verbatim
        Chen, J. H., Pei, A. G., Chen, P., & Hu, Z. Q. (2021).
        Study on gyroscopic effect of floating offshore wind turbines.
        China Ocean Engineering, 35(2), 201-214.
    \endverbatim

    Example specification in dynamicMeshDict:
    \verbatim
     restraints
     {
      turbineGyroscopic
      {
        sixDoFRigidBodyMotionRestraint gyroscopicMoment;
        rotationAxis      (1 0 0); // Initial rotation axis (in body frame)
        angularSpeed      1; // Angular speed around the specified axis, in [rad/s]            
        inertiaMoment     (1 0 0 0 1 0 0 0 1); //Initial moment of inertia tensor of the rotating body (in body frame)
      }
     }
    \endverbatim

SourceFiles
    gyroscopicMoment.C

Author
    Pere Frontera Pericàs, TUDelft
\*---------------------------------------------------------------------------*/

#ifndef gyroscopicMoment_H
#define gyroscopicMoment_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class gyroscopicMoment Declaration
\*---------------------------------------------------------------------------*/

class gyroscopicMoment
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data

        //- Initial rotation axis in body frame
        vector axis_;

        //- Angular speed in [rad/s]
        scalar speed_;

        //- Angular velocity (computed from axis_ and speed_)
        vector omega_;

        //- Initial moment of inertia of the rotating body in body frame
        tensor inertia_;

public:

    //- Runtime type information
    TypeName("gyroscopicMoment");


    // Constructors

        //- Construct from components
        gyroscopicMoment
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new gyroscopicMoment(*this)
            );
        }


    //- Destructor
    virtual ~gyroscopicMoment();


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

        //- Compute rotation matrix from axis-angle notation
        Foam::tensor rotAxis2Matrix(const vector& axis,const scalar& angle) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
