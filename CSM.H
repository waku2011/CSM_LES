/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

Class
    Foam::LESModels::CSM

Group
    grpLESTurbulence

Description
    The Coherent Structure model by H.Kobayashi(2005).

    Reference:
    \verbatim
        Kobayashi, H. The subgrid-scale models based on coherent structures
        for rotating homogeneous turbulence and turbulent channel flow, 
        Phys. Fluids, 17(2005),045104
    \endverbatim

    The default model coefficients are
    \verbatim
        CSMCoeffs
        {
            C2                  1/22;
        }
    \endverbatim

See also
    Foam::LESModels::Smagorinsky
    Foam::LESModels::WALE

SourceFiles
    CSM.C

\*---------------------------------------------------------------------------*/

#ifndef CSM_H
#define CSM_H

#include "LESModel.H"
#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class CSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class CSM
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        CSM(const CSM&);
        void operator=(const CSM&);

protected:

    // Protected data

        dimensionedScalar Ce_;
        dimensionedScalar C2_;

    // Protected Member Functions

        //- Return the deviatoric symmetric part of the square of the given
        //  velocity gradient field
        tmp<volSymmTensorField> Sd(const volTensorField& gradU) const;
  
        //- Return Choherent Structure Function, Fcs 
        tmp<volScalarField> Q(const volTensorField& gradU) const;
        tmp<volScalarField> E(const volTensorField& gradU) const;
 
        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const volTensorField& gradU) const;

        //- Update the SGS eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("CSM");


    // Constructors

        //- Construct from components
        CSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~CSM()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(this->U_));
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
