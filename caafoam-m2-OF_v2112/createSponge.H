    const vectorField& Cf = mesh.C();
 
    const scalar pi = constant::mathematical::pi;
//    dimensionedScalar r1 ("r1", dimensionSet(0, 1, 0, 0, 0),  300 );

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


   /*forAll(DT.boundaryField(),patchi)
   {
      if (isA<wallFvPatch>(mesh.boundary()[patchi]))
      {
        DT.boundaryField()[patchi] == 0;
      }
      else
      {
        DT.boundaryField()[patchi] == DT.boundaryField()[patchi].patchInternalField();
      }
   }*/

   forAll(blendFactor_.boundaryField(),patchi)
   {
      if ( mesh.boundary()[patchi].type() == "patch" || mesh.boundary()[patchi].type() == "symmetry"  )
      //if ( mesh.boundary()[patchi].type() == "wall"  )
      {
          blendFactor_.boundaryFieldRef()[patchi] = yy;
      }    
   }

   /*label PatchID = mesh.boundaryMesh().findPatchID("in");
   blendFactor_.boundaryField()[PatchID] = yy;  

   label PatchID2 = mesh.boundaryMesh().findPatchID("out");
   blendFactor_.boundaryField()[PatchID2] = yy;  

   label PatchID3 = mesh.boundaryMesh().findPatchID("simm_up");
   blendFactor_.boundaryField()[PatchID3] = yy;  

   label PatchID4 = mesh.boundaryMesh().findPatchID("simm_dw");
   blendFactor_.boundaryField()[PatchID4] = yy;*/
  
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
    //dimensionedScalar("one", dimensionSet(1, -3, 0, 0, 0), 1.4)
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
       scalar d = mag(Cf[celli]);

       if ( alpha[celli] > Ls)
       {
            blendFactor_[celli] = 0.0;
       } 
       else
       { 
            scalar xx = sigma0*Foam::pow((Ls- alpha[celli])/Ls , 2.0) ;
            //scalar xx = 0.34*Foam::pow((20.0- alpha[celli])/20.0 , 2.0) ;
            blendFactor_[celli] = xx;
       }

    }

 



