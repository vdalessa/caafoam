
    IOdictionary spongeDict
    (
        IOobject
        (
            "spongeDict",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );

        scalar Ls(readScalar(spongeDict.lookup("Ls"))); 
        scalar f(readScalar(spongeDict.lookup("f"))); 
        scalar a(readScalar(spongeDict.lookup("a"))); 
        scalar M(readScalar(spongeDict.lookup("M"))); 
        scalar LAD(readScalar(spongeDict.lookup("LAD"))); 

        scalar gamma(readScalar(spongeDict.lookup("gamma"))); 
        scalar pinf(readScalar(spongeDict.lookup("pinf"))); 
        scalar Tinf(readScalar(spongeDict.lookup("Tinf"))); 
        
        scalar Uinf_x(readScalar(spongeDict.lookup("Uinf_x"))); 
        scalar Uinf_y(readScalar(spongeDict.lookup("Uinf_y"))); 
        scalar Uinf_z(readScalar(spongeDict.lookup("Uinf_z"))); 


        Info <<  endl;
        
        Info<< "Ls = " << Ls << endl;
        Info<< "f = " << f << endl;
        Info<< "a = " << a << endl;
        Info<< "M = " << M << endl;
        Info<< "LAD = " << LAD << endl;
        
        Info <<  endl;
        Info<< "gamma = " << gamma << endl;
        Info<< "pinf  = " << pinf << endl;
        Info<< "Tinf  = " << Tinf << endl;

        Info <<  endl;
    //}
        scalar RT = a*a/gamma ;  
        scalar rho0 = pinf/RT ; 

        scalar Eref0 = (a*a/(gamma*(gamma-1.0))) + 0.5*( pow(Uinf_x,2) + pow(Uinf_y,2) + pow(Uinf_z,2) ) ; 

        Info<< "rho0  = " << rho0 << endl;
        Info<< "Eref0  = " << Eref0 << endl;
        Info<< "Uinf_x  = " << Uinf_x << endl;
        Info<< "Uinf_y  = " << Uinf_y << endl;
        Info<< "Uinf_z  = " << Uinf_z << endl;
        Info <<  endl;
     
        scalar sigma0 = (3.0*(1.0 - M*M))/(4.34294481903252e-01*Ls);
        Info<< "sigma0  = " << sigma0 << endl;
        Info <<  endl;
