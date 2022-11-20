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

#include "turbineAL.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(turbineAL, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        turbineAL,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::turbineAL::turbineAL
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict)
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::turbineAL::~turbineAL()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::turbineAL::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    //- Reference to time database
    const Time& time = motion.time();

    //- Create and initialise IO dictionary for turbine loads
    //- Only done at first time-step
    createTurbineDict(time);

    //- Access the dictionary
    dictionary loadDict = dictionary();

    // If turbine load IOdictionary exists
    if(time.foundObject<IOdictionary>("turbineSixDoFLoads"))
    {
        //- Access turbine load dict from the 'turbineSixDoFLoads' IOdictionary
        loadDict = time.lookupObject<IOdictionary>("turbineSixDoFLoads");
    }

    // Get load data from dictionary
    loadDict.lookup("refPoint") >> restraintPosition;
    loadDict.lookup("force") >> restraintForce;
    loadDict.lookup("moment") >> restraintMoment;

    if (motion.report())
    {
        Info<< " turbine force application point: " << restraintPosition << endl
            << " turbine force vector: " << restraintForce << endl
            << " turbine moment vector: " << restraintMoment << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::turbineAL::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::turbineAL::write
(
    Ostream& os
) const
{

}

void Foam::sixDoFRigidBodyMotionRestraints::turbineAL::createTurbineDict(const Time& time) const
{
    // Create dictionary if it has not been created before
    if(!time.foundObject<IOdictionary>("turbineSixDoFLoads"))
    {
        dictionary loadsDict;
        time.store
        (   
            new IOdictionary
            (   
                IOobject
                (   
                    "turbineSixDoFLoads",
                    time.timeName(),
                    time,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                loadsDict
            )
        );

        Info << "Rigid body IOdictionary 'turbineSixDoFLoads' created" << endl;

        initialiseTurbineDict(time);
    }
}

void Foam::sixDoFRigidBodyMotionRestraints::turbineAL::initialiseTurbineDict(const Time& time) const
{
    if(time.foundObject<IOdictionary>("turbineSixDoFLoads"))
    {
        // Open and write
        const dictionary& loadsDict = 
            time.lookupObject<IOdictionary>("turbineSixDoFLoads");

        // Initialise with zero entries
        dictionary updateDbDictionary = &loadsDict;
        updateDbDictionary.set("force", vector::zero);
        updateDbDictionary.set("moment", vector::zero);
        updateDbDictionary.set("refPoint", vector::zero);

        // Needed to update the IOdictionary
        const_cast<dictionary& > (loadsDict) = updateDbDictionary;
        
        Info << "Initialising 'turbineSixDoFLoads' IOdictionary" << endl;
    }
}

// ************************************************************************* //
