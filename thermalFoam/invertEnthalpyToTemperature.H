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
    invertEnthalpyToTemperature.H - H → T Inversion for Enthalpy-Primary Method

    PURPOSE
    -------
    Given enthalpy H, recover temperature T using the relationship:
        H = C_sens * T + f_l(T) * L * θ

    This is the critical step that allows the enthalpy equation to be solved
    directly while maintaining temperature for diffusion and boundary conditions.

    SUPPORTED SFCC MODELS
    ---------------------
    - "linear"   : Analytical inversion (piecewise linear)
    - "powerLaw" : Newton iteration required
    - "gaussian" : Newton iteration required

    THEORY
    ------
    For enthalpy-porosity method, we solve:
        ∂H/∂t = ∇·(K∇T)

    The H→T inversion handles the phase change singularity properly:
    - In frozen zone: T = H / C_sens
    - In mushy zone: T ≈ T_melt (clamped by latent heat)
    - In thawed zone: T = (H - L*θ) / C_sens

    The key insight is that during melting, H increases continuously while
    T stays nearly constant. This allows enthalpy diffusion to propagate
    heat into the interior even when temperature gradients are small.

    REFERENCES
    ----------
    - Voller & Swaminathan (1991) - Fixed grid enthalpy method
    - Alexiades & Solomon (1993) - Mathematical modeling of melting

\*---------------------------------------------------------------------------*/

#ifndef invertEnthalpyToTemperature_H
#define invertEnthalpyToTemperature_H

#include "fvCFD.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Invert H to T for LINEAR SFCC model (analytical)
inline scalar invertHtoT_linear
(
    const scalar H,
    const scalar C_sens,
    const scalar theta,
    const scalar L,
    const scalar T_sol,
    const scalar T_liq,
    const scalar thetar_frac
)
{
    // Handle dry cells
    if (theta < SMALL || C_sens < SMALL)
    {
        return H / max(C_sens, SMALL);
    }

    const scalar dT_mushy = max(T_liq - T_sol, SMALL);

    // Compute enthalpy bounds for phase transitions
    // H_sol: enthalpy at T = T_sol with f_l = thetar_frac
    // H_liq: enthalpy at T = T_liq with f_l = 1
    const scalar H_sol = C_sens * T_sol + thetar_frac * L * theta;
    const scalar H_liq = C_sens * T_liq + L * theta;

    if (H <= H_sol)
    {
        // Fully frozen: f_l = thetar_frac (constant)
        // H = C_sens * T + thetar_frac * L * theta
        return (H - thetar_frac * L * theta) / C_sens;
    }
    else if (H >= H_liq)
    {
        // Fully thawed: f_l = 1 (constant)
        // H = C_sens * T + L * theta
        return (H - L * theta) / C_sens;
    }
    else
    {
        // Mushy zone: f_l = thetar_frac + (1 - thetar_frac) * (T - T_sol) / dT_mushy
        // H = C_sens * T + [thetar_frac + (1 - thetar_frac) * (T - T_sol) / dT_mushy] * L * theta
        // H = C_sens * T + thetar_frac * L * theta + (1 - thetar_frac) * L * theta * (T - T_sol) / dT_mushy
        // H = T * [C_sens + (1 - thetar_frac) * L * theta / dT_mushy]
        //     + thetar_frac * L * theta - (1 - thetar_frac) * L * theta * T_sol / dT_mushy
        //
        // Solving for T:
        const scalar C_eff_mushy = C_sens + (1.0 - thetar_frac) * L * theta / dT_mushy;
        const scalar H_offset = thetar_frac * L * theta
                              - (1.0 - thetar_frac) * L * theta * T_sol / dT_mushy;

        return (H - H_offset) / C_eff_mushy;
    }
}


//- Invert H to T for POWER-LAW SFCC model (Newton iteration)
inline scalar invertHtoT_powerLaw
(
    const scalar H,
    const scalar C_sens,
    const scalar theta,
    const scalar L,
    const scalar T_f,
    const scalar T_ref,
    const scalar b_exp,
    const scalar thetar_frac,
    const scalar T_guess,
    const int maxIter = 20,
    const scalar tol = 1e-6
)
{
    // Handle dry cells
    if (theta < SMALL || C_sens < SMALL)
    {
        return H / max(C_sens, SMALL);
    }

    // Enthalpy at full thaw (T >= T_f, f_l = 1)
    const scalar H_thaw = C_sens * T_f + L * theta;

    if (H >= H_thaw)
    {
        // Above melting: f_l = 1
        return (H - L * theta) / C_sens;
    }

    // Newton iteration to find T such that H(T) = H_target
    // H(T) = C_sens * T + f_l(T) * L * theta
    // where f_l(T) = (1 - thetar_frac) * (1 - u(T)) + thetar_frac
    // and u(T) = ((T_f - T) / T_ref)^b for T < T_f

    scalar T = T_guess;

    for (int iter = 0; iter < maxIter; iter++)
    {
        scalar f_l, df_dT;
        if (T >= T_f)
        {
            f_l = 1.0; df_dT = 0.0;
        }
        else
        {
            const scalar depression = max(T_f - T, SMALL);
            const scalar b_val = max(mag(b_exp), SMALL);
            const scalar arg = depression / T_ref;
            const scalar r = Foam::pow(arg, b_val);
            
            // f_l = thetar + (1-thetar)/(1+r)
            f_l = thetar_frac + (1.0 - thetar_frac) / (1.0 + r);
            
            // df/dT = (1-thetar) * (b/Tr) * (dT/Tr)^(b-1) / (1+r)^2
            // Since dT/dT = -1 and df/dr is negative, df/dT is positive
            df_dT = (1.0 - thetar_frac) * (b_val / T_ref) 
                  * Foam::pow(arg, b_val - 1.0) / Foam::sqr(1.0 + r);
            
            df_dT = min(max(df_dT, scalar(0.0)), scalar(10.0));
        }

        // Compute H(T) and dH/dT
        const scalar H_calc = C_sens * T + f_l * L * theta;
        const scalar dH_dT = C_sens + df_dT * L * theta;

        // Newton update
        const scalar residual = H_calc - H;

        if (mag(residual) < tol * max(mag(H), scalar(1.0)))
        {
            break;  // Converged
        }

        if (mag(dH_dT) > SMALL)
        {
            T -= residual / dH_dT;
        }
        else
        {
            // Fallback: use sensible-only estimate
            T = H / C_sens;
        }
    }

    return T;
}


//- Invert H to T for GAUSSIAN SFCC model (Newton iteration)
inline scalar invertHtoT_gaussian
(
    const scalar H,
    const scalar C_sens,
    const scalar theta,
    const scalar L,
    const scalar T_f,
    const scalar W,
    const scalar thetar_frac,
    const scalar T_guess,
    const int maxIter = 20,
    const scalar tol = 1e-6
)
{
    // Handle dry cells
    if (theta < SMALL || C_sens < SMALL)
    {
        return H / max(C_sens, SMALL);
    }

    // Enthalpy at full thaw (T >= T_f, f_l = 1)
    const scalar H_thaw = C_sens * T_f + L * theta;

    if (H >= H_thaw)
    {
        // Above melting: f_l = 1
        return (H - L * theta) / C_sens;
    }

    // Newton iteration
    // f_l(T) = (1 - thetar_frac) * exp(-((T - T_f)/W)^2) + thetar_frac for T < T_f

    scalar T = T_guess;

    for (int iter = 0; iter < maxIter; iter++)
    {
        // Compute f_l and df_l/dT at current T
        scalar f_l, df_dT;

        if (T >= T_f)
        {
            f_l = 1.0;
            df_dT = 0.0;
        }
        else
        {
            const scalar arg = (T - T_f) / W;
            const scalar expTerm = Foam::exp(-arg * arg);

            f_l = (1.0 - thetar_frac) * expTerm + thetar_frac;

            // df_l/dT = (1 - thetar_frac) * exp(-arg^2) * (-2*arg/W)
            //         = -2 * (T - T_f) / W^2 * (1 - thetar_frac) * exp(-arg^2)
            df_dT = -2.0 * arg / W * (1.0 - thetar_frac) * expTerm;
        }

        // Compute H(T) and dH/dT
        const scalar H_calc = C_sens * T + f_l * L * theta;
        const scalar dH_dT = C_sens + df_dT * L * theta;

        // Newton update
        const scalar residual = H_calc - H;

        if (mag(residual) < tol * max(mag(H), scalar(1.0)))
        {
            break;  // Converged
        }

        if (mag(dH_dT) > SMALL)
        {
            T -= residual / dH_dT;
        }
        else
        {
            T = H / C_sens;
        }
    }

    return T;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
