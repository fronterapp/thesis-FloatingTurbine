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
    write1DFvMesh

Description

Author
    Niels Gjoel Jacobsen, Deltares

\*---------------------------------------------------------------------------*/

#include "write1DFvMesh.H"

#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

write1DFvMesh::write1DFvMesh
(
    const fvMesh& mesh,
    const word name
)
:
    mesh_(mesh),

    oneDFvMeshPtr_(NULL),

    polyPatches_(1),

    meshName_(name)
{
}


write1DFvMesh::~write1DFvMesh()
{
    if (oneDFvMeshPtr_ != NULL)
    {
        // Deleting pointers
        delete(oneDFvMeshPtr_);
        delete(polyPatches_[0]);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void write1DFvMesh::populatePoints
(
    const List<pointField>& lpp,
    const labelList& nPoints,
    pointField& pp
) const
{
    pp.setSize(Foam::sum(nPoints), point::zero);

    label offset = 0;

    forAll (lpp, listi)
    {
        const pointField& lp = lpp[listi];

        forAll (lp, pointi)
        {
            pp[pointi + offset] = lp[pointi];
        }

        offset += nPoints[listi];
    }
}


faceList write1DFvMesh::allFaces
(
    const label& nFaces,
    const label& nCells
) const
{
    faceList ff(nFaces);

    forAll (ff, facei)
    {
        face& f = ff[facei];
        f.setSize(4);
    }

    label C = 0;
    // Internal faces
    for (int i = 1; i < nCells; i++)
    {
        ff[C][0] = i*4;
        ff[C][1] = i*4 + 1;
        ff[C][2] = i*4 + 2;
        ff[C++][3] = i*4 + 3;
    }

    // First the end faces
    ff[C][0] = 0; ff[C][1] = 3; ff[C][2] = 2; ff[C++][3] = 1;

    ff[C][0] = nCells*4 + 0; ff[C][1] = nCells*4 + 1;
    ff[C][2] = nCells*4 + 2; ff[C++][3] = nCells*4 + 3;

    // Remaining faces - there is four faces around per cell
    for (int i = 0; i < nCells; i++)
    {
        // First face on circumference
           ff[C][0] = i*4;
           ff[C][1] = i*4 + 1;
           ff[C][2] = (i + 1)*4 + 1;
           ff[C++][3] = (i + 1)*4;

        // Second face on circumference
           ff[C][0] = i*4 + 1;
           ff[C][1] = i*4 + 2;
           ff[C][2] = (i + 1)*4 + 2;
           ff[C++][3] = (i + 1)*4 + 1;

        // Third face on circumference
           ff[C][0] = i*4 + 2;
           ff[C][1] = i*4 + 3;
           ff[C][2] = (i + 1)*4 + 3;
           ff[C++][3] = (i + 1)*4 + 2;

        // Fourth face on circumference
           ff[C][0] = i*4 + 3;
           ff[C][1] = i*4;
           ff[C][2] = (i + 1)*4;
           ff[C++][3] = (i + 1)*4 + 3;
    }

    return ff;
}


cellList write1DFvMesh::allCells(const label& nCells) const
{
    cellList cc(nCells);

    forAll (cc, celli)
    {
        cc[celli].setSize(6, -1);
    }

    label C = 0;
    // First cell
    cc[C][0] = 0;
    cc[C][1] = nCells - 1;
    cc[C][2] = nCells + 1;
    cc[C][3] = nCells + 2;
    cc[C][4] = nCells + 3;
    cc[C++][5] = nCells + 4;

    // Middle cells
    for ( int i = 1; i < nCells - 1; i++)
    {
        cc[C][0] = i - 1;
        cc[C][1] = i;
        cc[C][2] = nCells + i*4 + 1;
        cc[C][3] = nCells + i*4 + 2;
        cc[C][4] = nCells + i*4 + 3;
        cc[C++][5] = nCells + i*4 + 4;
    }

    // Last cell
    cc[C][0] = nCells - 2;
    cc[C][1] = nCells;
    cc[C][2] = 5*(nCells - 1) + 2;
    cc[C][3] = 5*(nCells - 1) + 3;
    cc[C][4] = 5*(nCells - 1) + 4;
    cc[C++][5] = 5*(nCells - 1) + 5;

    return cc;
}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void write1DFvMesh::updateMesh
(
      const List<pointField>& lpp,
       const labelList& nPoints,
       const labelList& nFaces,
       const labelList& nBndFaces,
       const labelList& nCells
)
{
    // Nothing to be done, if the number of cells is zero
    if (Foam::sum(nCells) == 0)
    {
        return;
    }

    // If the point is NULL, construct the mesh from scratch. Otherwise,
    // simply move the points
    if (oneDFvMeshPtr_ == NULL)
    {
        // Create the pointField
        pointField pp;
        populatePoints(lpp, nPoints, pp);

        // Create the faceList
        faceList faces;
        faces.setSize(Foam::sum(nFaces));

        label offsetInt = 0;
        label offsetBnd = Foam::sum(nFaces) - Foam::sum(nBndFaces);
        label offset = 0;

        forAll (lpp, listi)
        {
            // Construct the faces for a 1D mesh
              faceList ff = allFaces(nFaces[listi], nCells[listi]);

              // Loop over internal faces
              label intFaces = nFaces[listi] - nBndFaces[listi];

              for (int i = 0; i < intFaces; i++)
              {
                  face& f = faces[offsetInt + i];
                  f.setSize(4);
                  f = ff[i];

                  forAll (f, pointi)
                  {
                      f[pointi] += offset;
                  }
              }

              offsetInt += intFaces;

              // Loop over boundary faces
              for (int i = 0; i < nBndFaces[listi]; i++)
              {
                  face& f = faces[offsetBnd + i];
                  f.setSize(4);
                  f = ff[i + intFaces];

                  forAll (f, pointi)
                  {
                      f[pointi] += offset;
                  }
              }

              // Indent offsets
            offsetBnd += nBndFaces[listi];
            offset += nPoints[listi];
        }

        // Create the cellList in the memory
        cellList cells;
        cells.setSize(Foam::sum(nCells));

        // Set their size to 6
        forAll (cells, celli)
        {
            cell& c(cells[celli]);
            c.setSize(6, -1);
        }

        // Reset the offsets
        offsetInt = 0;
        offsetBnd = Foam::sum(nFaces) - Foam::sum(nBndFaces);
        offset = 0; // Used to offset the cell indices

        forAll (lpp, listi)
        {
            cellList cc = allCells(nCells[listi]);

            label intFaces = nFaces[listi] - nBndFaces[listi];

            forAll (cc, celli)
            {
                cell& C = cells[celli + offset];
                C = cc[celli];

                forAll (C, facei)
                {
                    if (C[facei] < intFaces)
                    {
                        C[facei] += offsetInt;
                    }
                    else
                    {
                        C[facei] += offsetBnd - intFaces;
                    }
                }
            }

            offsetInt += intFaces;
            offsetBnd += nBndFaces[listi];
            offset += nCells[listi];
        }


        // Create the fvMesh
        oneDFvMeshPtr_ = new fvMesh
        (
            IOobject
            (
                meshName_,
                mesh_.time().constant(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            std::move(pp),
            std::move(faces),
            std::move(cells)
        );

        // Add the boundary conditions by hand. Trivial, since the ordering and
        // type of boundaries are known apriori (chosen only cyclic in theta)
        polyBoundaryMesh bm(oneDFvMeshPtr_->thisDb(), *oneDFvMeshPtr_, 1);

        // Create a list of polyPatches
        polyPatches_[0] = new polyPatch
        (
            "boundary",
            Foam::sum(nBndFaces),
            Foam::sum(nFaces) - Foam::sum(nBndFaces),
            0,
            bm,
            "patch"
        );
    }
    else
    {
        // Populate the point fields
        pointField pp;

        populatePoints(lpp, nPoints, pp);

        // Move (update) the points
        oneDFvMeshPtr_->movePoints(pp);
    }
}


void write1DFvMesh::write(bool force) const
{
    if
    (
        oneDFvMeshPtr_ != NULL &&
        (mesh_.time().outputTime() || mesh_.time().timeIndex() < 1 || force)
    )
    {
        // Control the location, where the permeable mesh is written
        if (mesh_.time().timeIndex() > 0)
        {
            Info << "Time Name" << endl;
            oneDFvMeshPtr_->setInstance(mesh_.time().timeName());
        }
        else
        {
            Info << "Constant" << endl;
            oneDFvMeshPtr_->setInstance(mesh_.time().constant());
        }

        oneDFvMeshPtr_->write();
    }
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
