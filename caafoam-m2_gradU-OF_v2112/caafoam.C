////  =========                 |
//  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
//   \\    /   O peration     |
//    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
//     \\/     M anipulation  |
//-------------------------------------------------------------------------------
/*License
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

Application
    caafoam-m2 

Description
    Density based solver for Compressible Navier--Stokes equation.
    Low--storage fourth--order explicit Runge--Kutta scheme is used for time
    integration. 
    The convective terms are discretized using Pirozzoli's scheme(JCP 2010),
    density based splitting. The scheme allows conservation of the kinetic
    energy in the inviscid incompressible limit.
    midPoint interpolation must be selected in fvScheme dictionary. 
    Viscous terms(Laplacians only) are evaluated directly, computing
    the face normal gradients.
    Sponge--layer type mon reflective boundary treatment is  employed 
    to avoid numerical spurious reflections at far--field  boundaries.
    The OpenFOAM labrary for turbulence models is included.
 
    Authors:
    Valerio D'Alessandro - v.dalessandro@univpm.it
    Matteo Falone - m.falone@pm.univpm.it 
    Reference:
    V. D'Alessandro, M. Falone, R. Ricci.  Direct computation of aeroacoustic fields 
    in laminar flows: solver development and assessment of wall temperature effects 
    on radiated sound around bluff bodies. Comput. & Fluids (2020)

    Acknowledgements:  
    Authors thank Davide Modesti for sharing his solver, rhoEnergyFoam, which was
    was used as starting point for caafoam-m2.
    rhoEnergyFoam reference:
    D. Modesti and S. Pirozzoli A high-fidelity solver for turbulent compressible flows on unstructured meshes. Comput. & Fluids (2017)
*/
//
//\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
//#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include <fstream>      // std::ofstream

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// AUSM useful functions
  float m1 (float rm,float sgn ) // Return M_1
   {
   float r ; 
   r = 0.5*(rm + sgn*mag(rm)) ;
   return r ; 
   }
///
  float m2 (float rm,float sgn ) // Return M_2
   {
   float r ; 
   r = sgn*0.25*(rm + sgn)*(rm + sgn) ;
   return r ; 
   }
///
  float p5 (float rm,float sgn ,float alpha) // Return p_5
   {
   float r ; 
   if (abs(rm)<1.)
   {
    r = m2(rm,sgn)*( (sgn*2-rm) - sgn*16*alpha*rm*m2(rm,-sgn) ) ;
   }
   else
   {
    r = 1./rm*m1(rm,sgn) ;
   }
   return r ; 
   }
///
  float m4 (float rm,float sgn ,float beta) // Return M_4
   {
   float r ; 
   if (abs(rm)<1)
   {
    r = m2(rm,sgn)*( 1 - sgn*16*beta*m2(rm,-sgn) ) ;
   }
   else
   {
    r = m1(rm,sgn) ;
   }
   return r ; 
   }
// Main
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
 
    #include "addCheckCaseOptions.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "variables.H"
    #include "readSponge.H"
  	//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;
    Info<< "Start Timing = " << runTime.clockTimeIncrement() << " s"
        << nl << endl;

    rhoOld  = rho; 
    rhoUOld = rhoU; 
    rhoEOld = rhoE; 
//
    wrhoOld  = wrho; 
    wrhoUOld = wrhoU; 
    wrhoEOld = wrhoE; 
    #include "createSponge.H"

    while (runTime.loop()) //Start time loop
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;


//     Saving quantities at preavious time step
//     RK Time step

       for (int cycle =0; cycle < RK4values.size(); cycle++)
       {
//     Speed of sound and Mach number
        c = Foam::sqrt(thermo.Cp()/thermo.Cv()/psi);
        Mach = U/c ;
//      Interpolated quantities at cell faces     
        surfaceScalarField rhoave = fvc:: interpolate(rho) ;
        surfaceVectorField Uave   = fvc:: interpolate(U)   ;
//      Flux at the intercell
        phi    = fvc:: interpolate(U)    & mesh.Sf() ;
        phit   = fvc:: interpolate(rhoU) & mesh.Sf() ; //flux in turbulent model
        phit.setOriented(true);
//      Enthalpy
        H    = (rhoE + p)/rho ;
//      Enthalpy at the intercell
        surfaceScalarField Have = fvc:: interpolate(H) ;
//      Pressure at the intercell
        surfaceScalarField pave = fvc::interpolate(p)  ;       
        #include "sensor.H" //Ducros sensor 
        if(pressArtDiff)
        {
         #include "AUSM.H"    // AUSM+up dissipation on pressure term, add dissipation on pave
        }
//
//      Evaluate viscous terms
//
//      Divergence of the velocity             
        surfaceVectorField duV =  -2./3.*fvc::interpolate(fvc::div(U))*mesh.Sf();
        duV.setOriented(true);                                


	volTensorField     gU     = fvc::grad(U) ;
        surfaceTensorField gU_ave = fvc::interpolate( gU.T()  ) ;
        surfaceVectorField gU_tp  =  gU_ave & (mesh.Sf())  ;
//
        volScalarField muEff(turbulence->muEff());
        surfaceScalarField  muave  = fvc::interpolate(muEff);//mu at cell faces
//
        volScalarField k("k", thermo.Cp()*muEff/0.71);//thermal diffusivity
//
        surfaceScalarField kave=fvc::interpolate(k);//k at cell faces. alphaEff=muEff/Prt
        //momentum viscous flux
        surfaceVectorField momVisFlux = muave*(fvc::snGrad(U)*mesh.magSf() + duV + gU_tp);
        momVisFlux.setOriented(true);
        //energy viscous flux
        surfaceScalarField heatFlux =  kave*fvc::snGrad(T)*mesh.magSf();
        surfaceScalarField visWork  =  momVisFlux & Uave;
        enVisFlux = heatFlux + visWork ;
//
        // Total fluxes, Eulerian + viscous 
        surfaceScalarField rhoFlux   = -rhoave*phi ;   
        rhoFlux.setOriented(false);                                 ;
        momFlux                      = -rhoave*Uave*phi - pave*mesh.Sf() + momVisFlux  ;
        enFlux                       = -rhoave*Have*phi + enVisFlux                    ;
//
        if(convArtDiff)
        {
         #include "AUSM_conv.H"    // AUSM+up dissipation on convective terms, add dissipation on rhoFlux,momFlux,enFlux
        }
//
// 
        volScalarField rhoFl = fvc::div(rhoFlux) + blendFactor_*(rhoRef - rho)  ;
        volVectorField momFl = fvc::div(momFlux) + blendFactor_*(rhoRef*URef - rho*U);
        volScalarField enFl  = fvc::div(enFlux)  + blendFactor_*(rhoRef*ERef - rhoE);
        // RK sub-step
        wrho = RK4values[cycle]*wrhoOld + runTime.deltaT()*(rhoFl);    
        rho  = rhoOld + wrho*RK4values2[cycle];

        wrhoU = RK4values[cycle]*wrhoUOld + runTime.deltaT()*(momFl);  
        rhoU  = rhoUOld +  wrhoU * RK4values2[cycle];

        //Update primitive variables and boundary conditions
        U.ref() =
                    rhoU()
                   /rho();

        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef()         = rho.boundaryField()*U.boundaryField();

        
        wrhoE = RK4values[cycle]*wrhoEOld + runTime.deltaT()*(enFl);  
        rhoE  = rhoEOld +  wrhoE * RK4values2[cycle];

        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        //Thermodinamic library
        thermo.correct();
        rhoE.boundaryFieldRef() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );
//
        p.ref() =
                    rho()
                   /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() = psi.boundaryField()*p.boundaryField(); //psi=1/(R*T)
        runTime.write();
        turbulence->correct(); //turbulence model
    
        rhoOld = rho;
        rhoEOld = rhoE;
        rhoUOld = rhoU;

        wrhoOld = wrho;
        wrhoEOld = wrhoE;
        wrhoUOld = wrhoU;
       
       }//end of RK time integration
//
//
//
          Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
 
    
    }
    runTime.write();
    Info<< "Start Timing = " << runTime.clockTimeIncrement() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
