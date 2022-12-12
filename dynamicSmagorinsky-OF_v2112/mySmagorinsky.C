/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "mySmagorinsky.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> mySmagorinsky<BasicTurbulenceModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    volSymmTensorField D(symm(gradU));

    volScalarField a(this->Ce_/this->delta());
    volScalarField b((2.0/3.0)*tr(D));
    volScalarField c(2*Ck_*this->delta()*(dev(D) && D));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a))
        )
    );
}

template<class BasicTurbulenceModel>
void mySmagorinsky<BasicTurbulenceModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    // The SGS viscosity is bounded so that nuEff cannot become negative.
    // This allows backscatter (small scales transferring energy to large
    // scales) by allowing nut() to become -ve by allowing nuEff to be zero
    // when nut() goes -ve.     
       volSymmTensorField S(symm(gradU));          // Traced Strain Rate Sij
       volScalarField magS(sqrt(2.0*(S && S)));    // Norm |S|= sqrt(2.Sij.Sij)
    //
       this->nut_ = max(cD(S)*sqr(this->delta())*magS,-1.0*this->nu());
       this->nut_.correctBoundaryConditions();
       fv::options::New(this->mesh_).correct(this->nut_);
    
     BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void mySmagorinsky<BasicTurbulenceModel>::correctNut()
{
    /*volScalarField k(this->k(fvc::grad(this->U_)));

    this->nut_ = Ck_*this->delta()*sqrt(k);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
   
    BasicTurbulenceModel::correctNut();*/

     correctNut(fvc::grad(this->U_));
}

template<class BasicTurbulenceModel>
volScalarField mySmagorinsky<BasicTurbulenceModel>::cD(const volSymmTensorField& S) const
{
    const rhoField& rho = this->rho_;
    const volVectorField& U = this->U_;

    volScalarField magS(sqrt(2.0*(S && S)));

    volScalarField rhof(filter_(rho));          // rhof = <rho>
    volVectorField Uf(filter_(rho*U)/rhof);     // <U> = <rho.U>/<rho>
    volSymmTensorField Sf(symm(fvc::grad(Uf))); // <Sij> = 0.5(grad(<U>)+gradT(<U>))
    volScalarField magSf(sqrt(2.0*(Sf && Sf))); // Norm |<Sij>| = sqrt(2.<Sij><Sij>)

    // Leonard Term - Lij = <rho.Ui.Uj> - <rho><Ui><Uj>
    //     // ------------------------------------------------
    volSymmTensorField Lij = filter_(rho*sqr(U)) - sqr(filter_(rho*U))/rhof;
    //

       volSymmTensorField Bij
       (
        -2.0*sqr(2.0*this->delta())*rhof*magSf*dev(Sf)
       );
       volSymmTensorField Aij
       (
         -2.0*sqr(this->delta())*rho*magS*dev(S)
       );
      
       volSymmTensorField Mij = Bij - filter_(Aij);
       // LijMij = |Lij.Mij|
       //     // ------------------
       volScalarField LijMij = fvc::average(Lij && Mij);

       //MklMkl = |Mij.Mij|
       // ------------------ 
       volScalarField MklMkl = fvc::average(magSqr(Mij));

       // Cs = |Lij.Mij|/|Mij.Mij|
       // ------------------------
       // Cs = LijMij/MklMkl;
       
       volScalarField Cs = Cf_*0.0;
       forAll(this->mesh_.cells(),celli)
       {
            if (MklMkl[celli] != 0.0)
                Cs[celli] = LijMij[celli]/MklMkl[celli];
            else
                Cs[celli] = 0.0;
       }

       return Cs;

}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mySmagorinsky<BasicTurbulenceModel>::mySmagorinsky
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

    /*Ck_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    ),*/

       Cf_ // MA 19.JAN.2020
    (
        IOobject
        (
            IOobject::groupName("dynCf", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar ("dynCf",dimless, scalar(0.0))
    ),

 
    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_())
   

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mySmagorinsky<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mySmagorinsky<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mySmagorinsky<BasicTurbulenceModel>::omega() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));
    volScalarField epsilon(this->Ce_*k*sqrt(k)/this->delta());

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        epsilon/(0.09*k)
    );
}


template<class BasicTurbulenceModel>
void mySmagorinsky<BasicTurbulenceModel>::correct()
{
//    LESeddyViscosity<BasicTurbulenceModel>::correct();
//    correctNut();

    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();
    volSymmTensorField S(symm(gradU));                  // Traced Strain Rate Sij
    volScalarField magS(sqrt(2.0*(S && S)));    // Norm |S|= sqrt(2.Sij.Sij)

    // Calc SGS Eddy Viscosity
      correctNut(gradU);
    //
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
