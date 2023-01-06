
    volScalarField blendFactor_
    (
    IOobject
    (
        "blendFactor_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 0, -1, 0, 0), 0.0)
   );
 
//-------------------------------------------------------------------- 
   scalar yy = sigma0;


   forAll(blendFactor_.boundaryField(),patchi)
   {
      if ( mesh.boundary()[patchi].type() == "patch" || mesh.boundary()[patchi].type() == "symmetry"  )
      {
          blendFactor_.boundaryFieldRef()[patchi] = yy;
      }    
   }

  
//-------------------------------------------------------------------- 
//-------------------------------------------------------------------- 
  volScalarField rhoRef
    (
    IOobject
    (
        "rhoRef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(1, -3, 0, 0, 0), rho0)
   );

    volVectorField URef
    (
    IOobject
    (
        "URef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("one", dimensionSet(0, 1, -1, 0, 0), vector (Uinf_x, Uinf_y, Uinf_z) )
   );

    volScalarField ERef
    (
    IOobject
    (
        "ERef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 2, -2, 0, 0), Eref0)
   );


    forAll(mesh.C(), celli)
    {

       if ( alpha[celli] > Ls)
       {
            blendFactor_[celli] = 0.0;
       } 
       else
       { 
            scalar xx = sigma0*Foam::pow((Ls- alpha[celli])/Ls , 2.0) ;
            blendFactor_[celli] = xx;
       }

    }
 



