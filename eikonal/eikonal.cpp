/*---------------------------------------------------------------------------*\
Application
    eikonal 

Description
    Description
    Numerical solver for Eikonal equation which is solved to obtain the minimum 
    distance from specified boundaries. 
    The solver is developed within OpenFOAM 2.3.x.

    For correspondence:
    Valerio D'Alessandro - v.dalessandro@univpm.it
    Matteo Falone - m.falone@pm.univpm.it 

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    dimensionedScalar one ( "one", dimensionSet(0, 0, 0, 0, 0, 0, 0), scalar(1.0 ) );
    dimensionedScalar one1 ( "one1", dimensionSet(0, 1, 0, 0, 0, 0, 0), scalar(1.0 ) );
    dimensionedScalar one2 ( "one2", dimensionSet(0, 2, 0, 0, 0, 0, 0), scalar(1.0 ) );
    dimensionedScalar one_t ( "one_t", dimensionSet(0, 0, 1, 0, 0, 0, 0), scalar(1.0 ) );
  
    dimensionedScalar Ls ( "Ls", dimensionSet(0, 1, 0, 0, 0, 0, 0), scalar(0.01 ) );

    volVectorField gradU = fvc::grad(alpha);
    surfaceScalarField phi= (fvc::interpolate(gradU) & mesh.Sf());
//      
    while (runTime.loop())
    {
        Info<< "Iteration = " << runTime.timeName() << nl << endl;

            fvScalarMatrix alphaEqn
            (
                   one_t*fvm::ddt(alpha) +  one2*fvc::div(phi ,alpha)  -one - one2*alpha*fvc::div(phi)   // -scalar(1.0)
            );

            alphaEqn.relax();
            alphaEqn.solve();
            gradU = fvc::grad(alpha);
            phi= (fvc::interpolate(gradU) & mesh.Sf());
        
            volScalarField aa = (mag( alpha - alpha.oldTime() )/mag(alpha) ) ;
            scalar res = gMax(aa) ; 
            Info<< "Res = " << res << nl << endl;

            runTime.write();

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


        alpha.write() ; 

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
