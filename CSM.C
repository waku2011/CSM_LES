/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "CSM.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> CSM<BasicTurbulenceModel>::Sd
(
    const volTensorField& gradU
) const
{
    // gradU = { [du/dx,du/dy,du/dz]
    //           [dv/dx,dv/dy,dv/dz]
    //           [dw/dx,dw/dy,dw/dz] } 
    // S = gradU & gradU is a tensor 
    // symm(S) is a symmetrical tensor of S
    //     ST = symm(T) = 0.5*(S + S.T())
    // dev(ST) is a ST - 1/3*(tr(ST))I 
    //     tr(ST) = ST11 + ST22 + ST33
    
    return dev(symm(gradU & gradU));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> CSM<BasicTurbulenceModel>::Q
(
    const volTensorField& gradU
) const
{
  // other way of calculation Q
  volSymmTensorField S(symm(gradU));  // symmetric part of tensor
  volTensorField W(skew(gradU));      // anti-symmetric part
  volScalarField SS(S && S);
  volScalarField WW(W && W);
  return 0.5*(WW - SS);
  //return 0.5*( sqr(tr(gradU))-tr((gradU & gradU)) );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> CSM<BasicTurbulenceModel>::E
(
    const volTensorField& gradU
) const
{
  // other way of calculation Q
  volSymmTensorField S(symm(gradU));  // symmetric part of tensor
  volTensorField W(skew(gradU));      // anti-symmetric part
  volScalarField SS(S && S);
  volScalarField WW(W && W);
  return 0.5*(WW + SS);
  //return 0.5*(gradU && gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> CSM<BasicTurbulenceModel>::k
(
    const volTensorField& gradU
) const
{
    volScalarField magSqrSd(magSqr(Sd(gradU)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            (2.0*C2_/Ce_)*sqr(this->delta())*magSqr(dev(symm(gradU)))
        )
    );
}

template<class BasicTurbulenceModel>
void CSM<BasicTurbulenceModel>::correctNut()
{

    //- Return Choherent Structure Function, Fcs = Q/E
    volScalarField Q(this->Q(fvc::grad(this->U_)));
    volScalarField E(this->E(fvc::grad(this->U_)));

    volScalarField Fcs = min(max(-1.0, Q/max(E,dimensionedScalar("VSMALL",E.dimensions(),VSMALL))), 1.0);
    
    // Info << "max/min of Fcs (|Fcs|<1);" << endl;
    // Info << " Qmax " << max(Q) << endl;
    // Info << " Qmin " << min(Q) << endl; 
    // Info << " Emax " << max(E) << endl;
    // Info << " Emin " << min(E) << endl; 
    // Info << " Fcsmax " << max(Fcs) << endl;
    // Info << " Fcsmin " << min(Fcs) << endl; 

    // Ccsm
    volScalarField Ccsm = C2_*pow(mag(Fcs),1.5)*(1.0-Fcs);  

    //Info << " Ccsmmax " << max(Ccsm) << endl;
    //Info << " Ccsmmin " << min(Ccsm) << endl; 

    this->nut_ = Ccsm * sqr(this->delta()) * mag(dev(symm(fvc::grad(this->U_))));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
CSM<BasicTurbulenceModel>::CSM
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.0/22.0
        )
    )

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool CSM<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        C2_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> CSM<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}


template<class BasicTurbulenceModel>
void CSM<BasicTurbulenceModel>::correct()
{
    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
