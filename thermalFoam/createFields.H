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

/*--------------------------------*- C++ -*----------------------------------*\
|   thermalFoam - createFields.H                                              |
|   v2509-source-implicit - Field initialization for source-based EP method   |
|                                                                             |
|   Key features:                                                             |
|   - H field for diagnostics/output (T is primary variable)                  |
|   - C_sens is primary heat capacity (no latent spike)                       |
|   - liquidFraction updated from T in solver, not from SFCC                  |
|   - No static caching of parameters                                         |
\*---------------------------------------------------------------------------*/

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// ==========================================================================
// PHASE PROPERTIES
// ==========================================================================

const dimensionedScalar thetas("thetas", transportProperties);
const dimensionedScalar Cthsoil("Cthsoil", transportProperties);

const dimensionedScalar theta_org("theta_org", transportProperties);
const dimensionedScalar Cthorg("Cthorg", transportProperties);

const dimensionedScalar Cthw("Cthw", transportProperties);
const dimensionedScalar Cthice("Cthice", transportProperties);
const dimensionedScalar Cthair("Cthair", transportProperties);
const dimensionedScalar L("L", transportProperties);

// ==========================================================================
// THERMAL CONDUCTIVITY CONSTANTS (Johansen Model)
// ==========================================================================

Info<< "Reading thermal conductivity constants for Johansen model..." << endl;
const dimensionedScalar Kthw("Kthw", dimPower/dimLength/dimTemperature, transportProperties);
const dimensionedScalar Kthice("Kthice", dimPower/dimLength/dimTemperature, transportProperties);
const dimensionedScalar Kthsoil("Kthsoil", dimPower/dimLength/dimTemperature, transportProperties);
const dimensionedScalar Kthair("Kthair", dimPower/dimLength/dimTemperature, transportProperties);

Info<< "  Kthw   = " << Kthw.value() << " W/m/K" << endl;
Info<< "  Kthice = " << Kthice.value() << " W/m/K" << endl;
Info<< "  Kthsoil= " << Kthsoil.value() << " W/m/K" << endl;
Info<< "  Kthair = " << Kthair.value() << " W/m/K" << endl;

// ==========================================================================
// SFCC (Soil Freezing Characteristic Curve) PARAMETERS
//
// Supported SFCC models:
//   "linear"   - Linear mushy zone (default, original thermalFoam behavior)
//   "powerLaw" - Power-law: f_l = 1 - ((Tmelt - T)/Tref)^b  (McKenzie et al.)
//   "gaussian" - Gaussian:  f_l = exp(-((T - Tmelt)/W)^2)   (permaFoam style)
//
// ==========================================================================

const dimensionedScalar Tmelt("Tmelt", transportProperties);
const dimensionedScalar SFCC_Tref("SFCC_Tref", transportProperties);
const dimensionedScalar SFCC_b("SFCC_b", transportProperties);
const dimensionedScalar SFCC_S("SFCC_S", dimless, 0.01);
const dimensionedScalar Tminc("Tminc", transportProperties);
const dimensionedScalar Tmaxc("Tmaxc", transportProperties);

// SFCC model selection
const word SFCCModel = transportProperties.lookupOrDefault<word>("SFCCModel", "linear");

// Gaussian SFCC parameter: W [K] - controls freezing curve width
const scalar SFCC_W = transportProperties.lookupOrDefault<scalar>("SFCC_W", 1.0);

// Residual unfrozen water fraction [-] (water that doesn't freeze)
const scalar thetar_frac = transportProperties.lookupOrDefault<scalar>("thetar_frac", 0.0);

// SFCC smoothing width [K] - used in legacy apparent heat capacity mode
const scalar SFCC_smoothWidth = transportProperties.lookupOrDefault<scalar>("SFCC_smoothWidth", 0.5);

Info<< "\n=== SFCC Model Configuration ===" << nl
    << "  Model type     : " << SFCCModel << nl;
if (SFCCModel == "powerLaw")
{
    Info<< "  T_ref          : " << SFCC_Tref.value() << " K" << nl
        << "  b (exponent)   : " << SFCC_b.value() << nl;
}
else if (SFCCModel == "gaussian")
{
    Info<< "  W (width)      : " << SFCC_W << " K" << nl;
}
else if (SFCCModel == "linear")
{
    Info<< "  (Uses T_solidus/T_liquidus from Picard dict)" << nl;
}
Info<< "  thetar_frac    : " << thetar_frac << nl
    << "  Tmelt          : " << Tmelt.value() << " K" << nl
    << "================================" << endl;


// ==========================================================================
// THAW-BASED EROSION PARAMETERS
// ==========================================================================

#include "createErosionControls.H"

// ==========================================================================
// PRIMARY FIELDS
// ==========================================================================

Info<< "Reading field T (Temperature)\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field theta (Total Water Content)\n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// ==========================================================================
// ENTHALPY FIELD (for diagnostics and output)
//
// In source-based EP method, T is primary variable. H is computed from:
//   H = C_sens * T + f_l * L * θ
// ==========================================================================

Info<< "Initializing H (Enthalpy) field\n" << endl;
volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimVolume, 0.0)
);

// ==========================================================================
// THERMAL PROPERTY FIELDS
// ==========================================================================

Info<< "Initializing Kth (thermal conductivity via Johansen model)\n" << endl;
volScalarField Kth
(
    IOobject
    (
        "Kth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Kth", dimPower/dimLength/dimTemperature, 1.5)
);

volScalarField Cth
(
    IOobject
    (
        "Cth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", Cthsoil.dimensions(), 0.0)
);

volScalarField Cdg
(
    IOobject
    (
        "Cdg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Cdg", Cthsoil.dimensions(), 0.0)
);

// ==========================================================================
// SENSIBLE HEAT CAPACITY (primary for source-based EP method)
//
// C_sens = sum of component heat capacities weighted by volume fraction
// Does NOT include latent heat - that's handled by the source term
// ==========================================================================

volScalarField C_sens
(
    IOobject
    (
        "C_sens",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", Cthsoil.dimensions(), 0.0)
);

// Apparent heat capacity (for legacy AHC method compatibility)
volScalarField C_app
(
    IOobject
    (
        "C_app",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", Cthsoil.dimensions(), 0.0)
);

// Effective heat capacity
// - Source-based EP: CthEff = C_sens
// - Legacy AHC: CthEff = C_app (includes latent spike)
volScalarField CthEff
(
    IOobject
    (
        "CthEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", Cthsoil.dimensions(), 0.0)
);

// ==========================================================================
// PHASE FRACTION FIELDS
// ==========================================================================

volScalarField alphaWater
(
    IOobject
    (
        "alphaWater",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Liquid fraction - KEY FIELD for enthalpy-porosity method
// Updated from T in the Picard loop, not from SFCC
volScalarField liquidFraction
(
    IOobject
    (
        "liquidFraction",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Liquid water content
volScalarField thetal
(
    IOobject
    (
        "thetal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    theta
);

// Ice content
volScalarField thetag
(
    IOobject
    (
        "thetag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Ice fraction (complement of liquid fraction)
volScalarField alphaIce
(
    IOobject
    (
        "alphaIce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);

// ==========================================================================
// EROSION TRACKING FIELDS
// ==========================================================================

// Erosion mask: marks cells for removal (0 = not eroded, 1 = eroded)
volScalarField erosionMask
(
    IOobject
    (
        "erosionMask",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField originalPatchID
(
    IOobject
    (
        "originalPatchID",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("none", dimless, -1)
);

// ==========================================================================
// INITIALIZE originalPatchID FOR BOUNDARY TRACKING
// ==========================================================================

if (!originalPatchID.headerOk())
{
    Info<< "Initializing originalPatchID field" << endl;
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        // Check if this patch is in the erosionBoundaries list
        if (erosionBoundaries.found(pp.name()))
        {
            const labelList& faceCells = pp.faceCells();
            forAll(faceCells, faceI)
            {
                originalPatchID[faceCells[faceI]] = patchI;
            }
        }
    }
}

// ==========================================================================
// INITIALIZE LIQUID FRACTION FROM TEMPERATURE
//
// Used when starting fresh or when liquidFraction not in time directory.
// Uses the selected SFCC model (linear, powerLaw, or gaussian).
// ==========================================================================

if (!liquidFraction.headerOk())
{
    Info<< "Initializing liquidFraction from temperature using " << SFCCModel << " SFCC..." << endl;

    // Parameters for each SFCC model
    const scalar T_f_init = Tmelt.value();
    const scalar T_ref_init = max(SFCC_Tref.value(), SMALL);
    const scalar b_init = SFCC_b.value();
    const scalar W_init = max(SFCC_W, SMALL);

    // Linear model uses solidus/liquidus
    const scalar T_sol_init = transportProperties.lookupOrDefault<scalar>
    (
        "T_solidus",
        Tmelt.value() - 1.0
    );
    const scalar T_liq_init = transportProperties.lookupOrDefault<scalar>
    (
        "T_liquidus",
        Tmelt.value()
    );
    const scalar dT_mushy_init = max(T_liq_init - T_sol_init, SMALL);

    forAll(T, cellI)
    {
        const scalar T_cell = T[cellI];
        scalar fl = 0.0;

        if (SFCCModel == "linear")
        {
            if (T_cell <= T_sol_init)
            {
                fl = thetar_frac;
            }
            else if (T_cell >= T_liq_init)
            {
                fl = 1.0;
            }
            else
            {
                fl = thetar_frac + (1.0 - thetar_frac) * (T_cell - T_sol_init) / dT_mushy_init;
            }
        }
        else if (SFCCModel == "powerLaw")
        {
            // Power-law: f_l = (1 - thetar_frac) * (1 - u) + thetar_frac
            // where u = ((Tmelt - T)/Tref)^b
            if (T_cell >= T_f_init)
            {
                fl = 1.0;
            }
            else
            {
                const scalar depression = T_f_init - T_cell;
                const scalar T_ref_val = max(T_ref_init, SMALL);
                const scalar b_val = max(mag(b_init), SMALL); // Ensure positive magnitude
                
                // Standard retention curve: fl = thetar + (1-thetar)/(1 + (dT/Tr)^b)
                const scalar arg = depression / T_ref_val;
                const scalar r = Foam::pow(arg, b_val);
                fl = thetar_frac + (1.0 - thetar_frac) / (1.0 + r);
            }
        }
        else if (SFCCModel == "gaussian")
        {
            // Gaussian: f_l = (1 - thetar_frac) * exp(-((T - Tmelt)/W)^2) + thetar_frac
            if (T_cell >= T_f_init)
            {
                fl = 1.0;
            }
            else
            {
                const scalar arg = (T_cell - T_f_init) / W_init;
                fl = (1.0 - thetar_frac) * Foam::exp(-arg * arg) + thetar_frac;
            }
        }

        liquidFraction[cellI] = fl;

        // Update phase fractions
        thetal[cellI] = theta[cellI] * liquidFraction[cellI];
        thetag[cellI] = theta[cellI] * (1.0 - liquidFraction[cellI]);
    }

    alphaIce = scalar(1.0) - liquidFraction;

    Info<< "  liquidFraction range: [" << min(liquidFraction).value()
        << ", " << max(liquidFraction).value() << "]" << endl;
}

// ==========================================================================
// INITIALIZE ENTHALPY FROM TEMPERATURE AND LIQUID FRACTION
// (Used when starting fresh or when H not in time directory)
// ==========================================================================

if (!H.headerOk())
{
    Info<< "Initializing enthalpy H from temperature and liquid fraction..." << endl;
    
    // First compute sensible heat capacity
    forAll(C_sens, cellI)
    {
        C_sens[cellI] = 
              ((scalar(1.0) - thetas.value() - theta_org.value()) * Cthsoil.value())
            + (theta_org.value() * Cthorg.value())
            + (thetal[cellI] * Cthw.value())
            + (thetag[cellI] * Cthice.value())
            + (max(thetas.value() - theta[cellI], scalar(0.0)) * Cthair.value());
    }
    
    // Then compute enthalpy: H = C_sens * T + f_l * L * theta
    forAll(H, cellI)
    {
        H[cellI] = C_sens[cellI] * T[cellI] 
                 + liquidFraction[cellI] * L.value() * theta[cellI];
    }
    
    Info<< "  H range: [" << min(H).value() << ", " << max(H).value() << "] J/m³" << endl;
}

// ==========================================================================
// END OF FIELD INITIALIZATION
// ==========================================================================
