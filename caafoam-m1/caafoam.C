/*---------------------------------------------------------------------------*\
Application
    caafoam-m1

Description
    Description
    Density based solver for Compressible Navier--Stokes equation.
    Low--storage fourth--order explicit Runge--Kutta scheme is used for time
    integration. Central-upwind scheme of Kurganov-Noelle-Petrova is 
    used for convective terms. Sponge--layer type mon reflective boundary
    treatment is  employed to avoid numerical spurious reflections at far--field
    boundaries.
    The solver is developed within OpenFOAM 2.3.x.

    For correspondence:
    Valerio D'Alessandro - v.dalessandro@univpm.it
    Matteo Falone - m.falone@pm.univpm.it 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "turbulenceModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readThermophysicalProperties.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"
    #include "readSponge.H" 
 
    Info << mesh.geometricD() << endl;
    Info << mesh.solutionD() << endl;
    Info << mesh.nGeometricD() << endl;
    Info << mesh.nSolutionD() << endl;

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;



    while (runTime.run())
    {


//	Info << RK4values2 << endl;

	volScalarField rhoOld("rhoOld",rho);
	volScalarField rhoEOld("rhoEOld",rhoE);
	volVectorField rhoUOld("rhoUOld",rhoU);	

	volScalarField wrhoOld("wrhoOld",wrho);
	volScalarField wrhoEOld("wrhoEOld",wrhoE);
	volVectorField wrhoUOld("wrhoUOld",wrhoU);	

	#include "constructFluxes.H" 
	#include "createSponge.H" 

  
	for( int cycle = 0; cycle < RK4values.size(); cycle++)
	{
                sigmaU = 
                         ( (( (2.0/3.0)*muEff*fvc::div(U))*I + 2.0*muEff*symm(fvc::grad(U)))
                           & U
                         );

                rhokcycle = -fvc::div(phi) + blendFactor_*(rhoRef - rho) ;
                rhoUkcycle  = -fvc::div(phiUp) + fvc::laplacian(muEff, U) + fvc::div(tauMC) + blendFactor_*(rhoRef*URef - rho*U);
                rhoEkcycle =
                (
                 - fvc::div(phiEp)
                 + fvc::div(sigmaU)
                 + fvc::laplacian(k, T)
                 + blendFactor_*(rhoRef*ERef - rhoE)
                );
		
		// --- Solve density
                wrho = RK4values[cycle]*wrhoOld + runTime.deltaT()*rhokcycle ; 
		rho  = rhoOld  + wrho * RK4values2[cycle];

		// --- Solve momentum
                wrhoU = RK4values[cycle]*wrhoUOld + runTime.deltaT()*rhoUkcycle ; 
                rhoU  = rhoUOld + wrhoU * RK4values2[cycle] ;

		U.dimensionedInternalField() =
		    rhoU.dimensionedInternalField()
		   /rho.dimensionedInternalField();
		U.correctBoundaryConditions();
		rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

		// --- Solve energy
                wrhoE = RK4values[cycle]*wrhoEOld + runTime.deltaT()*rhoEkcycle ; 
                rhoE  = rhoEOld + wrhoE * RK4values2[cycle] ;

		e = rhoE/rho - 0.5*magSqr(U);
		e.correctBoundaryConditions();
		thermo.correct();

		rhoE.boundaryField() =
		    rho.boundaryField()*
		    (
		        e.boundaryField() + 0.5*magSqr(U.boundaryField())
		    );	

		p.dimensionedInternalField() =
		    rho.dimensionedInternalField()
		   /psi.dimensionedInternalField();
		p.correctBoundaryConditions();
		rho.boundaryField() = psi.boundaryField()*p.boundaryField();

		#include "fluxes.H" 

                /*sigmaU = 
                         ( (( (2.0/3.0)*muEff*fvc::div(U))*I + 2.0*muEff*symm(fvc::grad(U)))
                           & U
                         );*/
                  
           	k = thermo.Cp()*muEff/Pr;

                rhoOld = rho;
		rhoEOld = rhoE;
		rhoUOld = rhoU;	

                wrhoOld = wrho;
		wrhoEOld = wrhoE;
		wrhoUOld = wrhoU;	


	}// end loop

// finished calcuating fluxes constructing the 
        turbulence->correct();
        runTime.write();

	//#include "residuals.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
