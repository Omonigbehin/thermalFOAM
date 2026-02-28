/*---------------------------------------------------------------------------*\
  Copyright (C) 2026 Olorunfemi Omonigbehin, INRS-ETE
  This file is part of thermalFOAM
    
enthalpyUpdate.H - H→T Recovery for Enthalpy-Primary Method

    v2506-enthalpy-primary: ENTHALPY-PRIMARY FORMULATION

    PURPOSE
    -------
    This file provides the H→T inversion (temperature recovery from enthalpy)
    for use in the enthalpy-primary method. It is included here as an alternative solution

    In the main solver (thermalFoam.C), this functionality is implemented
    inline using the invertEnthalpyToTemperature.H header functions.

    THEORY
    ------
    The enthalpy-primary method solves:
        ∂H/∂t = ∇·(K∇T)

    After solving for H, we must recover T from the relationship:
        H = C_sens * T + f_l(T) * L * θ

    This is a nonlinear equation in T (because f_l depends on T via SFCC).
    The inversion is model-dependent:
    - Linear SFCC: Analytical piecewise solution
    - Power-law/Gaussian SFCC: Newton iteration required

    SUPPORTED SFCC MODELS
    ---------------------
    - "linear"   : H = C_sens*T for frozen, H = C_sens*T + L*θ for thawed
                   Mushy zone: linear interpolation
    - "powerLaw" : f_l = (1-θr)*(1 - ((Tm-T)/Tr)^b) + θr
    - "gaussian" : f_l = (1-θr)*exp(-((T-Tm)/W)^2) + θr

    REFERENCES
    ----------
    - Voller & Swaminathan (1991) - Fixed grid enthalpy method
    - Alexiades & Solomon (1993) - Mathematical modeling of melting

\*---------------------------------------------------------------------------*/

// ==========================================================================
// H → T RECOVERY (using inversion functions)
//
// This section recovers temperature from enthalpy using the SFCC model.
// For use in alternative workflows or post-processing.
// ==========================================================================

label nConverged_htinv = 0;
label nIterations_htinv = 0;

forAll(T, cellI)
{
    const scalar H_cell = H[cellI];
    const scalar C_sens_cell = C_sens[cellI];
    const scalar theta_cell = theta[cellI];
    const scalar T_current = T[cellI];  // Current value as initial guess

    scalar T_new = T_current;

    if (SFCCModel == "linear")
    {
        T_new = invertHtoT_linear
        (
            H_cell,
            C_sens_cell,
            theta_cell,
            L.value(),
            T_sol,
            T_liq,
            thetar_frac
        );
        nConverged_htinv++;
    }
    else if (SFCCModel == "powerLaw")
    {
        T_new = invertHtoT_powerLaw
        (
            H_cell,
            C_sens_cell,
            theta_cell,
            L.value(),
            T_f,
            T_ref,
            b_exp,
            thetar_frac,
            T_current,
            20,     // maxIter
            1e-6    // tol
        );
        nConverged_htinv++;
    }
    else if (SFCCModel == "gaussian")
    {
        T_new = invertHtoT_gaussian
        (
            H_cell,
            C_sens_cell,
            theta_cell,
            L.value(),
            T_f,
            W_gauss,
            thetar_frac,
            T_current,
            20,     // maxIter
            1e-6    // tol
        );
        nConverged_htinv++;
    }

    T[cellI] = T_new;
}

T.correctBoundaryConditions();

// ==========================================================================
// T → LIQUID FRACTION (SFCC)
//
// Update liquid fraction to be consistent with recovered temperature.
// ==========================================================================

label nSolid_eu = 0;
label nMushy_eu = 0;
label nLiquid_eu = 0;

forAll(liquidFraction, cellI)
{
    const scalar T_cell = T[cellI];
    scalar fl_new = 0.0;

    if (SFCCModel == "linear")
    {
        if (T_cell <= T_sol)
        {
            fl_new = thetar_frac;
            nSolid_eu++;
        }
        else if (T_cell >= T_liq)
        {
            fl_new = 1.0;
            nLiquid_eu++;
        }
        else
        {
            fl_new = thetar_frac + (1.0 - thetar_frac) * (T_cell - T_sol) / dT_mushy;
            nMushy_eu++;
        }
    }
    else if (SFCCModel == "powerLaw")
    {
        if (T_cell >= T_f)
        {
            fl_new = 1.0;
            nLiquid_eu++;
        }
        else
        {
            const scalar depression = max(T_f - T_cell, scalar(0.0));
            const scalar b_val = max(mag(b_exp), SMALL);
            const scalar r = Foam::pow(max(depression / T_ref, SMALL), b_val);
            
            // Correct retention form mapping
            fl_new = thetar_frac + (1.0 - thetar_frac) / (1.0 + r);

            // Updated threshold for reporting
            const scalar solidThreshold_eu = thetar_frac + 0.001;
            
            if (fl_new <= solidThreshold_eu)
                nSolid_eu++;
            else if (fl_new >= 0.999)
                nLiquid_eu++;
            else
                nMushy_eu++;
        }
    }
    else if (SFCCModel == "gaussian")
    {
        if (T_cell >= T_f)
        {
            fl_new = 1.0;
            nLiquid_eu++;
        }
        else
        {
            const scalar arg = (T_cell - T_f) / W_gauss;
            fl_new = (1.0 - thetar_frac) * Foam::exp(-arg * arg) + thetar_frac;

            if (fl_new < 0.01)
                nSolid_eu++;
            else if (fl_new > 0.99)
                nLiquid_eu++;
            else
                nMushy_eu++;
        }
    }

    // Apply under-relaxation if specified
    liquidFraction[cellI] = flRelaxation * fl_new
                          + (1.0 - flRelaxation) * liquidFraction[cellI];

    // Enforce bounds
    liquidFraction[cellI] = min(max(liquidFraction[cellI], scalar(0.0)), scalar(1.0));
}

liquidFraction.correctBoundaryConditions();

// ==========================================================================
// UPDATE alphaIce (complement of liquidFraction)
// ==========================================================================

alphaIce = scalar(1.0) - liquidFraction;
alphaIce.correctBoundaryConditions();

// ==========================================================================
// RE-COMPUTE H FOR CONSISTENCY
//
// After updating T and f_l, ensure H = C_sens*T + f_l*L*θ exactly.
// ==========================================================================

forAll(H, cellI)
{
    H[cellI] = C_sens[cellI] * T[cellI]
             + liquidFraction[cellI] * L.value() * theta[cellI];
}

H.correctBoundaryConditions();

// ==========================================================================
// REPORT
// ==========================================================================

reduce(nSolid_eu, sumOp<label>());
reduce(nMushy_eu, sumOp<label>());
reduce(nLiquid_eu, sumOp<label>());

Info<< "    H→T recovery: " << nConverged_htinv << " cells inverted" << nl
    << "    Phase distribution: solid=" << nSolid_eu
    << ", mushy=" << nMushy_eu
    << ", liquid=" << nLiquid_eu << nl
    << "    f_l range=[" << min(liquidFraction).value()
    << ", " << max(liquidFraction).value() << "]"
    << nl;

// ************************************************************************* //
