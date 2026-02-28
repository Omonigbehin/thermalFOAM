/*---------------------------------------------------------------------------*\
  Copyright (C) 2025 Olorunfemi Omonigbehin, INRS-ETE

  This file is part of thermalFOAM.

  thermalFOAM is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  thermalFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with thermalFOAM. If not, see <https://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
    updateConstitutive.H - Unified Physics Update for Implicit Enthalpy Method

    v2512-wave-saturation: WAVE-DRIVEN PORE SATURATION MODEL

    This file computes:
    1. Phase fractions (thetal, thetag) from liquidFraction
       WITH wave-driven saturation: as ice melts, freed pore space AND
       existing air voids fill with liquid water from wave action.
    2. Thermal conductivity (Johansen geometric mean model)
    3. Sensible heat capacity C_sens (no latent contribution)
    4. Effective heat capacity CthEff = C_sens + L·θ_orig·|∂f_l/∂T| (APPARENT)
       NOTE: Latent heat uses ORIGINAL theta (only original ice content melts)

    KEY DESIGN PRINCIPLES
    ---------------------
    - liquidFraction is the SOURCE OF TRUTH for phase state
    - thetal, thetag are DERIVED from liquidFraction (one-way dependency)
    - As melting progresses (f_l increases), wave-driven water fills:
      (a) pore space freed by melting ice
      (b) air-filled voids that were initially present
    - This is a simplified surrogate for Richards-equation water infiltration
      (cf. permaFoam, Orgogozo et al., 2023, Comput. Phys. Commun.)
    - Latent heat ONLY acts on the original ice content (theta_orig)
      Extra liquid water arriving from waves is already liquid — no phase change
    - Kth and C_sens use the EFFECTIVE (saturating) phase fractions

    PHYSICS
    -------
    Initial state:
      theta_orig = 0.198  (ice + residual liquid)
      porosity   = 0.448
      air voids  = 0.448 - 0.198 = 0.250

    As f_l goes from thetar_frac -> 1.0:
      Ice melts:      theta_ice = theta_orig * (1 - f_l)  [decreasing]
      Freed space:    theta_orig * f_l  [becomes liquid]
      Air fills with water from waves, proportional to melting progress:
        theta_wave = (f_l - thetar) / (1 - thetar) * (porosity - theta_orig)
      Total liquid:   theta_orig * f_l + theta_wave
      Total water:    theta_orig + theta_wave  (ice + liquid, increases toward porosity)
      Air remaining:  porosity - theta_orig - theta_wave  (decreases toward 0)

    CALLED FROM
    -----------
    1. At initialization (t=0) in thermalFoam.C
    2. Every Picard iteration in the thermal solution loop
    3. After mesh changes (erosion events)

    NO STATIC VARIABLES - safe for multiple includes per timestep.

\*---------------------------------------------------------------------------*/

// ========================================================================
// 1. PHASE FRACTIONS with WAVE-DRIVEN PORE SATURATION
//
// As ice melts, wave-driven water fills the freed pore space AND the
// pre-existing air voids. The saturation fraction ramps linearly with
// liquid fraction from 0 (fully frozen) to 1 (fully thawed = fully
// saturated by wave water).
//
// Latent heat (computed in Section 4) uses theta_orig — only the
// original ice undergoes phase change. The extra wave water is
// already liquid.
// ========================================================================

forAll(thetal, cellI)
{
    const scalar theta_orig = theta[cellI];    // original (frozen) water content
    const scalar porosity   = thetas.value();  // total pore space

    if (theta_orig < SMALL)
    {
        // Dry cell: no water content
        thetal[cellI] = 0.0;
        thetag[cellI] = 0.0;
        continue;
    }

    // Clamp liquidFraction for safety
    const scalar f_l = min(max(liquidFraction[cellI], scalar(0.0)), scalar(1.0));

    // -------------------------------------------------------------------
    // Wave-driven saturation fraction:
    // 0 when fully frozen (f_l = thetar_frac), 1 when fully thawed (f_l = 1)
    // -------------------------------------------------------------------
    scalar satFrac = 0.0;
    if (f_l > thetar_frac + SMALL)
    {
        satFrac = (f_l - thetar_frac) / max(1.0 - thetar_frac, SMALL);
        satFrac = min(max(satFrac, scalar(0.0)), scalar(1.0));
    }

    // Extra liquid water volume from wave infiltration into air voids
    const scalar airVoids = max(porosity - theta_orig, scalar(0.0));
    const scalar theta_wave = satFrac * airVoids;

    // Ice content: only the ORIGINAL water content contributes ice
    // (wave water arrives as liquid, never freezes in this context)
    const scalar phi_ice = theta_orig * (1.0 - f_l);

    // Liquid content: melted original ice + residual liquid + wave water
    const scalar phi_liq = theta_orig * f_l + theta_wave;

    thetal[cellI] = phi_liq;
    thetag[cellI] = phi_ice;
}

thetal.correctBoundaryConditions();
thetag.correctBoundaryConditions();

// Update alphaIce as complement of liquidFraction
alphaIce = scalar(1.0) - liquidFraction;
alphaIce.correctBoundaryConditions();


// ========================================================================
// 2. THERMAL CONDUCTIVITY (Johansen Geometric Mean Model)
//
// K_eff = K_soil^φ_solid * K_water^φ_liquid * K_ice^φ_ice * K_air^φ_air
//
// Now uses the EFFECTIVE phase fractions from Section 1, which include
// wave-driven water filling air voids. As melting progresses:
//   - φ_ice decreases (ice melts)
//   - φ_liq increases (melt water + wave water)
//   - φ_air decreases (filled by wave water)
//   - Kth increases significantly (water replaces insulating air)
// ========================================================================

forAll(Kth, cellI)
{
    // Solid fraction (mineral soil + organic matter)
    const scalar porosity = max(thetas.value(), scalar(0.0));
    scalar phi_solid = 1.0 - porosity;
    phi_solid = min(max(phi_solid, scalar(0.0)), scalar(1.0));

    // Phase volumes from computed fractions (now include wave water)
    scalar phi_liq = max(thetal[cellI], scalar(0.0));
    scalar phi_ice = max(thetag[cellI], scalar(0.0));

    // Enforce pore capacity: liquid + ice cannot exceed porosity
    scalar theta_tot = phi_liq + phi_ice;
    if (theta_tot > porosity && theta_tot > SMALL)
    {
        const scalar scale = porosity / theta_tot;
        phi_liq *= scale;
        phi_ice *= scale;
        theta_tot = porosity;
    }

    // Air-filled porosity (remaining pore space — decreases as waves fill voids)
    scalar phi_air = porosity - theta_tot;
    phi_air = min(max(phi_air, scalar(0.0)), porosity);

    // Renormalize to sum exactly 1
    scalar sumPhi = phi_solid + phi_liq + phi_ice + phi_air;
    if (sumPhi < SMALL) sumPhi = 1.0;

    phi_solid /= sumPhi;
    phi_liq   /= sumPhi;
    phi_ice   /= sumPhi;
    phi_air   /= sumPhi;

    // Guard against pow(0, phi) - use small positive values
    const scalar Ks = max(Kthsoil.value(), scalar(1e-12));
    const scalar Kw = max(Kthw.value(),    scalar(1e-12));
    const scalar Ki = max(Kthice.value(),  scalar(1e-12));
    const scalar Ka = max(Kthair.value(),  scalar(1e-12));

    // Geometric mean model (Johansen)
    Kth[cellI] =
        Foam::pow(Ks, phi_solid) *
        Foam::pow(Kw, phi_liq)   *
        Foam::pow(Ki, phi_ice)   *
        Foam::pow(Ka, phi_air);

    // Minimum conductivity floor
    Kth[cellI] = max(Kth[cellI], scalar(1e-4));
}

Kth.correctBoundaryConditions();


// ========================================================================
// 3. SENSIBLE HEAT CAPACITY
//
// C_sens = sum of component heat capacities weighted by volume fraction
//
// Uses EFFECTIVE phase fractions (with wave water), same as Kth.
// This is the PRIMARY heat capacity for source-based EP method.
// Latent heat is handled separately in Section 4 using theta_orig.
// ========================================================================

forAll(C_sens, cellI)
{
    // Mineral soil volume fraction (1 - porosity - organic)
    const scalar phi_mineral = max(scalar(1.0) - thetas.value() - theta_org.value(), scalar(0.0));

    // Use EFFECTIVE phase fractions (include wave water)
    const scalar phi_liq = thetal[cellI];
    const scalar phi_ice = thetag[cellI];

    // Air-filled porosity: what is left of pore space after wave filling
    const scalar phi_pore_filled = phi_liq + phi_ice;
    const scalar phi_air = max(thetas.value() - phi_pore_filled, scalar(0.0));

    // Weighted average of heat capacities
    C_sens[cellI] =
          phi_mineral * Cthsoil.value()           // Mineral soil
        + theta_org.value() * Cthorg.value()      // Organic matter
        + phi_liq * Cthw.value()                  // Liquid water (melt + wave)
        + phi_ice * Cthice.value()                // Ice (original only)
        + phi_air * Cthair.value();               // Air (decreasing)
}

C_sens.correctBoundaryConditions();


// ========================================================================
// 4. EFFECTIVE HEAT CAPACITY (Apparent Heat Capacity)
//
// CthEff = C_sens + L * theta_ORIG * |∂f_l/∂T|
//
// CRITICAL: Latent heat uses the ORIGINAL theta (from the theta field),
// NOT the effective theta. Only the ice that was originally present
// undergoes phase change. The extra wave-driven water is already liquid
// and contributes NO latent heat.
//
// This ensures:
//   - Correct latent heat barrier (only melting original ice)
//   - Correct Kth and C_sens (including wave water for transport)
//   - No artificial latent heat inflation from wave water
// ========================================================================

if (useEnthalpyPorosity)
{
    // Compute apparent heat capacity with latent contribution
    // CthEff = C_sens + L * theta_ORIG * |df_l/dT|

    const scalar T_f = Tmelt.value();
    const scalar T_ref = max(SFCC_Tref.value(), SMALL);
    const scalar b = SFCC_b.value();
    const scalar W = max(SFCC_W, SMALL);

    // For linear SFCC, need solidus/liquidus
    const scalar T_sol_const = transportProperties.lookupOrDefault<scalar>
    (
        "T_solidus",
        Tmelt.value() - 1.0
    );
    const scalar T_liq_const = transportProperties.lookupOrDefault<scalar>
    (
        "T_liquidus",
        Tmelt.value()
    );
    const scalar dT_mushy_const = max(T_liq_const - T_sol_const, SMALL);

    forAll(CthEff, cellI)
    {
        const scalar T_val = T[cellI];

        // ORIGINAL theta for latent heat — NOT the effective (wave-filled) value
        const scalar thetaOrig = theta[cellI];

        // Start with sensible capacity (which DOES use effective fractions)
        scalar C_eff_cell = C_sens[cellI];

        if (thetaOrig > SMALL)
        {
            scalar df_dT_mag = 0.0;

            if (SFCCModel == "linear")
            {
                if (T_val > T_sol_const && T_val < T_liq_const)
                {
                    df_dT_mag = (1.0 - thetar_frac) / dT_mushy_const;
                }
            }
            else if (SFCCModel == "powerLaw")
            {
                if (T_val < T_f)
                {
                    const scalar b_val = max(mag(b), SMALL);
                    const scalar arg = max((T_f - T_val) / T_ref, SMALL);
                    const scalar r = Foam::pow(arg, b_val);

                    df_dT_mag = (1.0 - thetar_frac) * (b_val / T_ref)
                              * Foam::pow(arg, b_val - 1.0) / Foam::sqr(1.0 + r);

                    df_dT_mag = min(df_dT_mag, scalar(10.0));
                }
            }
            else if (SFCCModel == "gaussian")
            {
                if (T_val < T_f)
                {
                    const scalar arg = (T_val - T_f) / W;
                    const scalar expTerm = Foam::exp(-arg * arg);
                    df_dT_mag = mag(-2.0 * arg / W * (1.0 - thetar_frac) * expTerm);
                }
            }

            // Latent contribution uses ORIGINAL theta only
            C_eff_cell += L.value() * thetaOrig * df_dT_mag;
        }

        CthEff[cellI] = C_eff_cell;
    }

    // Apply minimum floor
    CthEff = max
    (
        CthEff,
        dimensionedScalar("minC", C_sens.dimensions(), minCthEff)
    );
}
else
{
    // ====================================================================
    // LEGACY APPARENT HEAT CAPACITY METHOD
    //
    // Adds latent heat contribution via SFCC derivative:
    // C_app = C_sens + L * |∂θ_ice/∂T|
    //
    // Uses ORIGINAL theta for latent heat, same principle as above.
    // ====================================================================

    const scalar T_f = Tmelt.value();
    const scalar T_ref = max(SFCC_Tref.value(), SMALL);
    const scalar b = SFCC_b.value();
    const scalar S = max(SFCC_smoothWidth, SMALL);

    C_app = C_sens;

    forAll(C_app, cellI)
    {
        const scalar T_val = T[cellI];

        // ORIGINAL theta for latent heat
        const scalar thetaOrig = theta[cellI];

        if (thetaOrig < SMALL) continue;

        scalar dThetaIce_dT_mag = 0.0;

        if (T_val < T_f)
        {
            const scalar depression = T_f - T_val;

            if (depression >= S)
            {
                const scalar arg = max(depression / T_ref, SMALL);
                dThetaIce_dT_mag =
                    thetaOrig * mag(b / T_ref) * Foam::pow(arg, b - 1.0);
            }
            else
            {
                scalar u_at_S = Foam::pow(S / T_ref, b);
                u_at_S = min(max(u_at_S, scalar(0.0)), scalar(1.0));
                dThetaIce_dT_mag = thetaOrig * (1.0 - u_at_S) / S;
            }
        }

        C_app[cellI] += L.value() * dThetaIce_dT_mag;

        // Cap to prevent numerical issues
        C_app[cellI] = min(C_app[cellI], scalar(1e9));
    }

    CthEff = max
    (
        C_app,
        dimensionedScalar("minC", C_sens.dimensions(), minCthEff)
    );
}

// Store sensible capacity in Cth for output
Cth = C_sens;
Cth.correctBoundaryConditions();
CthEff.correctBoundaryConditions();

// ************************************************************************* //