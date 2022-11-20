/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::catenaryShape

Description


SourceFiles
    catenaryShape.C

Author
    Niels Gjoel Jacobsen, Deltares.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef catenaryShape_H
#define catenaryShape_H

#include "fvMesh.H"

namespace Foam
{

class catenaryShape
{
private:

    // Private member variables
        point pos0_;
        point pos1_;
        scalar lineLength_;
        scalar mass_;
        scalar maxCatenaryLineLength_;
        bool catenaryFormulation_;

        vector g_;
        scalar gMag_;

        scalar x0_, y0_, x1_, y1_, span_;

        scalar exceedStiffness_;
        scalar exceedLength_;
        scalar exceedForce_;

        scalar h_;
        scalar k_;

        scalar b_, X_, Y_;

    // Private member functions
        void shapeParameter();

        scalar shapeParameter(const scalar, const scalar, const scalar);

        //- Disallow default bitwise copy construct
        catenaryShape(const catenaryShape &);

        //- Disallow default bitwise assignment
        void operator=(const catenaryShape &);


public:

        catenaryShape
        (
            const point,
            const point,
            const scalar,
            const scalar,
            const vector
        );

        // Evaluate
//        void points(pointField&, const scalar) const;

        void centreLine(pointField&) const;

        vector H0() const;
        vector H1() const;

        vector R0() const;
        vector R1() const;

        bool isUShaped() const;

        inline scalar xTop() const {return X_;};


};

} // End namespace Foam

#endif
