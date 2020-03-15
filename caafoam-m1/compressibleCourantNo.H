
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

if (mesh.nInternalFaces())
{
    surfaceScalarField amaxSfbyDelta
    (
        mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
    );

    CoNum = max(amaxSfbyDelta/mesh.magSf()).value()*runTime.deltaTValue();

    meanCoNum =
        (sum(amaxSfbyDelta)/sum(mesh.magSf())).value()
       *runTime.deltaTValue();
}

Info<< "Mean and max Courant Numbers = "
    << meanCoNum << " " << CoNum << endl;

// ************************************************************************* //
