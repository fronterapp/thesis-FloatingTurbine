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

\*---------------------------------------------------------------------------*/

#include "gyroscopicMoment.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(gyroscopicMoment, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        gyroscopicMoment,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::gyroscopicMoment::gyroscopicMoment
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict),
    axis_(),
    speed_(),
    omega_(),
    inertia_()
{
    read(sDoFRBMRDict);

    //- Make sure that axis_ is normalised
    scalar magnitude = mag(axis_);
    if (magnitude>VSMALL)
    { 
        axis_ /= magnitude;
    }
    else
    {
        axis_ = vector(1, 0, 0); 
        Info<< " gyroscopicMoment restraint: rotationAxis is a null vector. " << endl
            << " New rotationAxis: " << axis_ << endl
            << endl;
    }

    //- Use axis_ and speed_ to compute the angular velocity vector
    omega_ = speed_*axis_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::gyroscopicMoment::~gyroscopicMoment()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::gyroscopicMoment::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{     
    //- Reference to time database
    const Time& time = motion.time();

    //- Angle rotated since initial time in [rad]
    scalar angle = speed_ * (time.value()-time.startTime().value()); // theta = w * (t-t0)

    //- Compute rotation matrix based on the rotation axis_ and angle
    tensor rotMatrix = rotAxis2Matrix(axis_, angle);

    //- Rotate the inertia tensor accordingly
    tensor inertiaRotated = rotMatrix & (inertia_ & rotMatrix.T());

    //- Compute gyroscopic load in global frame
    //- gyro_local = - omegaRigidBody_local x (I_local w_local)
    //- gyro_global = Q gyro_local = omegaRigidBody_global x (Q (I_local w_local))
    vector gyroMoment = -motion.omega() ^ (motion.orientation() & (inertiaRotated & omega_));

    // Applied force: none
    restraintForce = vector::zero;
    restraintPosition = vector::zero;

    // Applied moment
    restraintMoment = gyroMoment;  

    if (motion.report())
    {
        Info<< " Body angular velocity (body frame): " << omega_ << endl
            << " Body rotation angle (body frame): " << rotMatrix << endl
            << " Body inertia (body frame): " << inertiaRotated << endl
            << " Gyroscopic moment (global frame): " << restraintMoment << endl
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::gyroscopicMoment::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("rotationAxis", axis_);
    sDoFRBMRCoeffs_.readEntry("angularSpeed", speed_);
    sDoFRBMRCoeffs_.readEntry("inertiaMoment", inertia_);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::gyroscopicMoment::write
(
    Ostream& os
) const
{
    os.writeEntry("rotationAxis", axis_);
    os.writeEntry("angularSpeed", speed_);
    os.writeEntry("inertiaMoment", inertia_);
}

Foam::tensor Foam::sixDoFRigidBodyMotionRestraints::gyroscopicMoment::rotAxis2Matrix
(
    const vector& axis,
    const scalar& angle
) const
{
    // Declare and define the rotation matrix (from SOWFA)
    tensor RM;
    RM.xx() = Foam::sqr(axis.x())
            + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle);
    RM.xz() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y())
            + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z())
            + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    return RM;
}

// ************************************************************************* //
