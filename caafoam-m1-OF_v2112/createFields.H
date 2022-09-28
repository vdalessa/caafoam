Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();

volScalarField& p = thermo.p();
volScalarField& e = thermo.he();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

//bool inviscid(true);
//if (max(mu.internalField()) > 0.0)
//{
//    inviscid = false;
//}

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField wrho = rho;

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);
volVectorField wrhoU = rhoU;

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

volScalarField wrhoE = rhoE;

volScalarField alpha 
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);




surfaceScalarField phi("phi", mesh.Sf() & fvc::interpolate(rhoU));

// create the k values for the fields

scalarList RK4values(0);
scalarList RK4values2(0);

//RK4values.append(-1.0e9);
RK4values.append(0.0);
RK4values.append(-0.4178904745);
RK4values.append(-1.192151694643);
RK4values.append(-1.697784692471);
RK4values.append(-1.514183444257);

//RK4values2.append(-1.0e9);
RK4values2.append(0.1496590219993);
RK4values2.append(0.3792103129999);
RK4values2.append(0.8229550293869);
RK4values2.append(0.6994504559488);
RK4values2.append(0.1530572479681);

Info << RK4values2 << endl;
Info << RK4values << endl;

List<autoPtr<volScalarField> >  rhok(RK4values.size());
List<autoPtr<volVectorField> >  rhoUk(RK4values.size());
List<autoPtr<volScalarField> >  rhoEk(RK4values.size());

List<autoPtr<volScalarField> >  wrhok(RK4values.size());
List<autoPtr<volVectorField> >  wrhoUk(RK4values.size());
List<autoPtr<volScalarField> >  wrhoEk(RK4values.size());

volVectorField sigmaU
(
    IOobject
    (
        "sigmaU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
     mesh,
     dimensionedVector("zero", dimensionSet(1, 0, -3, 0, 0, 0, 0), vector(0.0, 0.0, 0.0))
);

 
volScalarField rhokcycle 
(
    IOobject
    (
        "rhokcycle",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
);
volVectorField rhoUkcycle 
(
    IOobject
    (
        "rhoUkcycle",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
     mesh,
     dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector(0.0, 0.0, 0.0))
);
volScalarField rhoEkcycle 
(
    IOobject
    (
        "rhoEkcycle",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      mesh,
      dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)
);



Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
