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
    Foam::sixDoFRigidBodyMotionRestraints::turbineAL

Description
    sixDoFRigidBodyMotionRestraints model.

    This model accounts for the loads from a turbine, 
    modelled using turbinesFOAM library. The 'coupleLoads' 
    flag in fvOptions must be set to true.

    The loads are stored in the 'turbineSixDoFLoads' IOdictionary, 
    which is created by this model. 

    Example specification in dynamicMeshDict:
    \verbatim
     restraints
     {
        turbineLoad
        {
            sixDoFRigidBodyMotionRestraint turbineAL; 
        }
     }
    \endverbatim

SourceFiles
    turbineAL.C

Author
    Pere Frontera Pericàs, TUDelft
\*---------------------------------------------------------------------------*/

#ifndef turbineAL_H
#define turbineAL_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class turbineAL Declaration
\*---------------------------------------------------------------------------*/

class turbineAL
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data
       
public:
    //- Runtime type information
    TypeName("turbineAL");

    // Public data

    // Constructors

        //- Construct from components
        turbineAL
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new turbineAL(*this)
            );
        }


    //- Destructor
    virtual ~turbineAL();


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

        //- Create the 'turbineSixDoFLoads' IOdictionary where turbinesFOAM will write turbine loads
        void createTurbineDict(const Time& time) const;

        //- Initialise 'turbineSixDoFLoads' IOdictionary
        void initialiseTurbineDict(const Time& time) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
