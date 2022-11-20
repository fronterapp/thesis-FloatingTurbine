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
    Foam::write1DFvMesh

Description
    Makes a 1D mesh along a centre line.


SourceFiles
    write1DFvMesh.C

Author
    Niels Gjoel Jacobsen, Deltares.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef write1DFvMesh_H
#define write1DFvMesh_H

#include "fvMesh.H"
#include "Time.H"

namespace Foam
{

class write1DFvMesh
{
private:

    // Private member variables
        const fvMesh& mesh_;

        fvMesh* oneDFvMeshPtr_;

        List<polyPatch*> polyPatches_;

        word meshName_;

    // Private member functions

        void populatePoints
        (
            const List<pointField>&,
            const labelList&,
            pointField&
        ) const;

        faceList allFaces(const label&, const label&) const;

        cellList allCells(const label&) const;

        //- Disallow default bitwise copy construct
        write1DFvMesh(const write1DFvMesh &);

        //- Disallow default bitwise assignment
        void operator=(const write1DFvMesh &);


public:

        write1DFvMesh
        (
            const fvMesh&,
            const word
        );

        ~write1DFvMesh();

        //- Create the 1D mesh from points field and cell information
        void updateMesh
        (
            const List<pointField>&,
            const labelList&,
            const labelList&,
            const labelList&,
            const labelList&
        );

        // Write
        void write(bool force = false) const;

        const fvMesh& oneDMesh() {return *oneDFvMeshPtr_;};


};

} // End namespace Foam

#endif
