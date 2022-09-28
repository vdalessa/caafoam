        // --- upwind interpolation of primitive fields on faces

        rho_pos = fvc::interpolate(rho, pos, "reconstruct(rho)");
        
        rho_neg =  fvc::interpolate(rho, neg, "reconstruct(rho)");

        rhoU_pos = fvc::interpolate(rhoU, pos, "reconstruct(U)");

        rhoU_neg = fvc::interpolate(rhoU, neg, "reconstruct(U)");


        rPsi = 1.0/psi;


        rPsi_pos = fvc::interpolate(rPsi, pos, "reconstruct(T)");

        rPsi_neg = fvc::interpolate(rPsi, neg, "reconstruct(T)");

        e_pos = fvc::interpolate(e, pos, "reconstruct(T)");

        e_neg = fvc::interpolate(e, neg, "reconstruct(T)");

        U_pos = rhoU_pos/rho_pos;
        U_neg = rhoU_neg/rho_neg;

        p_pos = rho_pos*rPsi_pos;
        p_neg = rho_neg*rPsi_neg;

        phiv_pos = U_pos & mesh.Sf();
        phiv_neg = U_neg & mesh.Sf();

        c = sqrt(thermo.Cp()/thermo.Cv()*rPsi);
        cSf_pos = fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf();

        cSf_neg = fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf();

        ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);

        am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);

        a_pos = ap/(ap - am);

        amaxSf = max(mag(am), mag(ap));

        aSf =  (am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        a_neg = 1.0 - a_pos;

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        aphiv_pos = (phiv_pos - aSf);
        aphiv_neg = (phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));





        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        phiUp =
		(
		    (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
		  + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
		);

        phiEp = 
		(
		    aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
		  + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
		  + aSf*p_pos - aSf*p_neg
		);

        volScalarField muEff(turbulence->muEff());
        tauMC = muEff*dev2(Foam::T(fvc::grad(U)));

