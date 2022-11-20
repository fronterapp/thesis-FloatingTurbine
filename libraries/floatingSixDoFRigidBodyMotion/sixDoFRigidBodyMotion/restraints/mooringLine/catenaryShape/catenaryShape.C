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
    catenaryShape

Description

Author
    Niels Gjoel Jacobsen, Deltares

\*---------------------------------------------------------------------------*/

#include "catenaryShape.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

catenaryShape::catenaryShape
(
    const point pos0,
    const point pos1,
    const scalar L,
    const scalar m,
    const vector gravity
)
:
    pos0_(pos0),
    pos1_(pos1),

    lineLength_(L),

    mass_(m),

    maxCatenaryLineLength_(0.999*lineLength_),

    catenaryFormulation_(true)
{
    g_ = gravity;
    gMag_ = Foam::mag(g_);
    g_ /= Foam::mag(g_);

    // Set the local coordinate system
    y0_ = (pos0_ & g_);
    y1_ = (pos1_ & g_);

    x0_ = 0.0;

    vector vec = pos1_ - pos0_;
    x1_ = Foam::mag(vec - (vec & g_)*g_);
    span_ = x1_;

    h_ = y1_ - y0_;

    if (maxCatenaryLineLength_ < Foam::sqrt(Foam::sqr(span_) + Foam::sqr(h_)))
    {
        catenaryFormulation_ = false;

        // Computing the shape parameter for two lengths of the line
        scalar span = Foam::sqrt(Foam::sqr(0.998*lineLength_) - Foam::sqr(h_));

        scalar k0 = shapeParameter(lineLength_, h_, span);

        span = Foam::sqrt(Foam::sqr(maxCatenaryLineLength_) - Foam::sqr(h_));
        scalar k1 = shapeParameter(lineLength_, h_, span);

        exceedForce_ = mass_*gMag_/k1;
        exceedLength_ = maxCatenaryLineLength_;
        exceedStiffness_ = (mass_*gMag_/k1 - mass_*gMag_/k0)/(0.001*lineLength_);

        b_ = -GREAT;
        X_ = -GREAT;
        Y_ = -GREAT;
    }
    else
    {
        shapeParameter();

        b_ = 2.0/k_*0.5*Foam::log((lineLength_ + h_)/(lineLength_ - h_));
        X_ = 0.5*(span_ + b_);
        Y_ = 1/k_*(Foam::cosh(k_*X_) - 1);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void catenaryShape::shapeParameter()
{
    k_ = shapeParameter(lineLength_, h_, span_);
}


scalar catenaryShape::shapeParameter
(
    const scalar lineLength,
    const scalar h,
    const scalar span
)
{
    scalar k = Foam::sqrt
        (
            24*(Foam::sqrt(Foam::sqr(lineLength) - Foam::sqr(h)) - span)
            /Foam::pow(span, 3.0)
        );

    while (true)
    {
        scalar f = k/2.0*sqrt(Foam::sqr(lineLength) - Foam::sqr(h))
            - Foam::sinh(span*k/2.0);

        if (Foam::mag(f) < 1.0e-5)
        {
            break;
        }

        scalar df = 1/2.0*sqrt(Foam::sqr(lineLength) - Foam::sqr(h))
            - span/2.0*Foam::cosh(span*k/2.0);

        k = k - f/df;
    }

    return k;
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void catenaryShape::centreLine(pointField& cl) const
{
    label N = cl.size();

    if (catenaryFormulation_)
    {
        scalarField x(N, 0.0), y(N, 0.0);
        scalar dx = span_/static_cast<scalar>(N - 1);

        forAll (x, i)
        {
            x[i] = i*dx;
            y[i] = 1/k_*(1 - cosh(k_*(x[i] - X_))) + Y_;
        }

        vector horzVec = pos1_ - pos0_;
        horzVec = horzVec - (horzVec & g_)*g_;
        horzVec /= Foam::mag(horzVec);

        cl = pos0_ + y*g_ + x*horzVec;
    }
    else
    {
        vector vec = pos1_ - pos0_;

        forAll (cl, i)
        {
            cl[i] = pos0_ + static_cast<scalar>(i)/static_cast<scalar>(N)*vec;
        }
    }
}


vector catenaryShape::H0() const
{
    vector vec = pos1_ - pos0_;
    vec /= Foam::mag(vec);

    if (catenaryFormulation_)
    {
        vec -= (vec & g_)*g_;

        vec *= mass_*gMag_/k_;
    }
    else
    {
        scalar dist = Foam::mag(pos1_ - pos0_);

        vec *= (exceedStiffness_*(dist - exceedLength_) + exceedForce_);
    }

    return vec;
}


vector catenaryShape::H1() const
{
    vector vec = pos1_ - pos0_;
    vec /= Foam::mag(vec);

    if (catenaryFormulation_)
    {
        vec -= (vec & g_)*g_;

        vec *= -mass_*gMag_/k_;
    }
    else
    {
        scalar dist = Foam::mag(pos1_ - pos0_);

        vec *= -(exceedStiffness_*(dist - exceedLength_) + exceedForce_);
    }

    return vec;
}


vector catenaryShape::R0() const
{
    if (catenaryFormulation_)
    {
        scalar H = mass_*gMag_/k_;
        scalar T = H*Foam::cosh(k_*(x0_ - X_));

        return g_*Foam::sqrt(Foam::sqr(T) - Foam::sqr(H));
    }
    else
    {
        return vector::zero;
    }
}


vector catenaryShape::R1() const
{
    if (catenaryFormulation_)
    {
        scalar H = mass_*gMag_/k_;
        scalar T = H*Foam::cosh(k_*(x1_ - X_));

        return g_*Foam::sqrt(Foam::sqr(T) - Foam::sqr(H));
    }
    else
    {
        return vector::zero;
    }
}


bool catenaryShape::isUShaped() const
{
    // The catenary is overstretched, i.e. it cannot be U-shaped.
    if (!catenaryFormulation_)
    {
        return false;
    }

    if (x0_ < X_ && X_ < x1_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
