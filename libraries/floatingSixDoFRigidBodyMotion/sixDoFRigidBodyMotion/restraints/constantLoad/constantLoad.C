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

#include "constantLoad.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(constantLoad, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        constantLoad,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::constantLoad::constantLoad
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict),
    applicationPt_(),
    force_(),
    torque_(),
    movePt_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::constantLoad::~constantLoad()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::constantLoad::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    if (movePt_)
    {
        // The force application point must move with the rigid body
        restraintPosition = motion.transform(applicationPt_); 
    }
    else
    {
        restraintPosition = applicationPt_;
    }

    // Applied force
    restraintForce = force_;

    // Applied torque
    restraintMoment = torque_;  

    if (motion.report())
    {
        Info<< " force application point: " << restraintPosition << endl
            << " force vector: " << restraintForce << endl
            << " torque vector: " << restraintMoment << endl
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::constantLoad::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("applicationPt", applicationPt_);
    sDoFRBMRCoeffs_.readEntry("force", force_);
    sDoFRBMRCoeffs_.readEntry("torque", torque_);
    sDoFRBMRCoeffs_.readEntry("movePt", movePt_);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::constantLoad::write
(
    Ostream& os
) const
{
    os.writeEntry("applicationPt", applicationPt_);
    os.writeEntry("force", force_);
    os.writeEntry("torque", torque_);
    os.writeEntry("movePt", movePt_);
}

// ************************************************************************* //
