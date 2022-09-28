		// --- upwind interpolation of primitive fields on faces

		surfaceScalarField rho_pos
		(
		    fvc::interpolate(rho, pos, "reconstruct(rho)")
		);
		surfaceScalarField rho_neg
		(
		    fvc::interpolate(rho, neg, "reconstruct(rho)")
		);

		surfaceVectorField rhoU_pos
		(
		    fvc::interpolate(rhoU, pos, "reconstruct(U)")
		);
		surfaceVectorField rhoU_neg
		(
		    fvc::interpolate(rhoU, neg, "reconstruct(U)")
		);

		volScalarField rPsi(1.0/psi);
		surfaceScalarField rPsi_pos
		(
		    fvc::interpolate(rPsi, pos, "reconstruct(T)")
		);
		surfaceScalarField rPsi_neg
		(
		    fvc::interpolate(rPsi, neg, "reconstruct(T)")
		);

		surfaceScalarField e_pos
		(
		    fvc::interpolate(e, pos, "reconstruct(T)")
		);
		surfaceScalarField e_neg
		(
		    fvc::interpolate(e, neg, "reconstruct(T)")
		);

		surfaceVectorField U_pos(rhoU_pos/rho_pos);
		surfaceVectorField U_neg(rhoU_neg/rho_neg);

		surfaceScalarField p_pos(rho_pos*rPsi_pos);
		surfaceScalarField p_neg(rho_neg*rPsi_neg);

		surfaceScalarField phiv_pos(U_pos & mesh.Sf());
                phiv_pos.setOriented(false);
		surfaceScalarField phiv_neg(U_neg & mesh.Sf());
                phiv_neg.setOriented(false);

		volScalarField c(sqrt(thermo.Cp()/thermo.Cv()*rPsi));
		surfaceScalarField cSf_pos
		(
		    fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf()
		);
		surfaceScalarField cSf_neg
		(
		    fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf()
		);

		surfaceScalarField ap
		(
		    max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
		);
		surfaceScalarField am
		(
		    min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
		);

		surfaceScalarField a_pos(ap/(ap - am));

		surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

		surfaceScalarField aSf(am*a_pos);

		if (fluxScheme == "Tadmor")
		{
		    aSf = -0.5*amaxSf;
		    a_pos = 0.5;
		}

		surfaceScalarField a_neg(1.0 - a_pos);

		phiv_pos *= a_pos;
		phiv_neg *= a_neg;

		surfaceScalarField aphiv_pos(phiv_pos - aSf);
		surfaceScalarField aphiv_neg(phiv_neg + aSf);

		// Reuse amaxSf for the maximum positive and negative fluxes
		// estimated by the central scheme
		amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

		#include "compressibleCourantNo.H"
		#include "readTimeControls.H"
		#include "setDeltaT.H"

		runTime++;

		Info<< "Time = " << runTime.timeName() << nl << endl;

		phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;
                surfaceVectorField phiU(aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg);
                phiU.setOriented(true);
       
                surfaceVectorField phiUp(phiU + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf());


		/*phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

		surfaceVectorField phiUp
		(
		    (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
		  + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
		);*/

		surfaceScalarField phiEp
		(
		    aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
		  + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
		  + aSf*p_pos - aSf*p_neg
		);

		volScalarField muEff(turbulence->muEff());
		volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

           	volScalarField k("k", thermo.Cp()*muEff/0.71);
           	//volScalarField k("k", thermo.Cp()*muEff/Pr);
