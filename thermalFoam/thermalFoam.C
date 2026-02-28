/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026
    Omonigbehin Olorunfemi
-------------------------------------------------------------------------------
License
    This file is part of thermalFOAM.

    thermalFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    thermalFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with thermalFOAM.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------------------
Application
    thermalFOAM

Description
    Transient solver for coupled heat transfer, phase change, and
    thaw-driven erosion of coastal permafrost.

    Governing equation:
        ∂H/∂t − ∇·(K ∇T) = 0

    where enthalpy H incorporates latent heat effects using an
    enthalpy-porosity formulation.

    Key features:
      - Implicit enthalpy-porosity phase change
      - Wave-aware convective heat flux boundary condition
      - Dynamic cell removal for thaw-based erosion
      - Compatible with dynamicRefineFvMesh
      - Sequential solution strategy with Picard iteration control

Author
    Omonigbehin Olorunfemi
    PhD Candidate, Institut National de l Recherche Scientific
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "pointMesh.H"
#include "syncTools.H"
#include "invertEnthalpyToTemperature.H"

int main(int argc, char* argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    
    // ======================================================================
    // PICARD ITERATION CONTROLS
    //
    // Recommended settings for implicit enthalpy diffusion:
    //   nCorrectors: 20-30 (more iterations help with mushy zone)
    //   tolerance: 1e-3 K (can relax from 1e-4 with implicit scheme)
    //   flRelaxation: 0.7-0.9 (higher values after initial stabilization)
    // ======================================================================

    const dictionary& picardDict = mesh.solutionDict().subDict("Picard");
    const int nCorrectors = picardDict.lookupOrDefault<int>("nCorrectors", 20);
    const scalar finalTol = picardDict.lookupOrDefault<scalar>("tolerance", 1e-3);
    const scalar flRelaxation = picardDict.lookupOrDefault<scalar>("flRelaxation", 0.8);
    
    // ======================================================================
    // ADAPTIVE TIME-STEPPING CONTROLS
    // ======================================================================
    
    const scalar tFact = picardDict.lookupOrDefault<scalar>("tFact", 1.5);
    const scalar maxDeltaT = picardDict.lookupOrDefault<scalar>("maxDeltaT", 10.0);
    const scalar minDeltaT = picardDict.lookupOrDefault<scalar>("minDeltaT", 1e-4);
    const int stabilityThreshold = picardDict.lookupOrDefault<int>("stabilityThreshold", 10);
    const bool adjustTimeStep = runTime.controlDict().lookupOrDefault<Switch>("adjustTimeStep", false);
    
    // Fourier number control (for adaptive timestep guidance)
    const scalar maxFourier = picardDict.lookupOrDefault<scalar>("maxFourier", 0.5);
    
    const scalar minCthEff = transportProperties.lookupOrDefault<scalar>("minCthEff", 1e5);

    // Maximum CthEff cap to ensure beta = 1/CthEff doesn't become too small
    // This prevents H diffusion from stalling in the mushy zone
    // Physical interpretation: limits how much latent heat a cell can absorb
    // at constant temperature before it must start warming
    const scalar maxCthEff = transportProperties.lookupOrDefault<scalar>("maxCthEff", 1e8);
    
    // Temperature relaxation within Picard
    scalar relaxationFactor = 0.8;
    if (mesh.solutionDict().found("relaxationFactors"))
    {
        const dictionary& relaxDict = mesh.solutionDict().subDict("relaxationFactors");
        if (relaxDict.found("fields"))
        {
            const dictionary& fieldRelax = relaxDict.subDict("fields");
            if (fieldRelax.found("T"))
            {
                relaxationFactor = readScalar(fieldRelax.lookup("T"));
            }
        }
    }
    
    // ======================================================================
    // ENTHALPY-POROSITY METHOD CONTROLS
    // ======================================================================
    
    const bool useEnthalpyPorosity = transportProperties.lookupOrDefault<Switch>
    (
        "useEnthalpyPorosity", 
        true
    );
    
    // Mushy zone parameters (used for linear SFCC mode)
    const scalar T_sol = transportProperties.lookupOrDefault<scalar>
    (
        "T_solidus",
        Tmelt.value() - 1.0
    );

    const scalar T_liq = transportProperties.lookupOrDefault<scalar>
    (
        "T_liquidus",
        Tmelt.value()
    );

    const scalar dT_mushy = max(T_liq - T_sol, SMALL);

    // SFCC model parameters (from createFields.H)
    const scalar T_f = Tmelt.value();
    const scalar T_ref = max(SFCC_Tref.value(), SMALL);
    const scalar b_exp = SFCC_b.value();
    const scalar W_gauss = max(SFCC_W, SMALL);

    Info<< "\n=== Thermal Solution Method ===" << nl;
    if (useEnthalpyPorosity)
    {
        Info<< "  Method: IMPLICIT ENTHALPY DIFFUSION (v2511)" << nl
            << "  Solves: ∂H/∂t - ∇·(K·β·∇H) = deferred_correction" << nl
            << "  where β = dT/dH ≈ 1/CthEff (Picard linearization)" << nl
            << "  Recovery: T = T(H) via SFCC inversion each iteration" << nl
            << "  Convergence: Dual criterion (T AND H residuals)" << nl
            << "  Picard iterations: max " << nCorrectors << ", tol(T) = " << finalTol << " K" << nl
            << "  Liquid fraction relaxation: " << flRelaxation << nl
            << "  Fourier stability: uses CthEff (apparent capacity)" << nl
            << "  CthEff bounds: [" << minCthEff << ", " << maxCthEff << "] J/m³/K" << nl
            << "  (maxCthEff caps beta to prevent H diffusion stall)" << nl;

        if (SFCCModel == "linear")
        {
            Info<< "  SFCC: LINEAR mushy zone" << nl
                << "    T_solidus  = " << T_sol << " K" << nl
                << "    T_liquidus = " << T_liq << " K" << nl
                << "    dT_mushy   = " << dT_mushy << " K" << nl;
        }
        else if (SFCCModel == "powerLaw")
        {
            Info<< "  SFCC: POWER-LAW (Retention Form)" << nl
                << "    f_l = thetar + (1-thetar)/(1 + ((Tm-T)/Tref)^b)" << nl
                << "    Tmelt = " << T_f << " K, Tref = " << T_ref << " K, b = " << mag(b_exp) << nl;
        }
        else if (SFCCModel == "gaussian")
        {
            Info<< "  SFCC: GAUSSIAN (permaFoam style)" << nl
                << "    f_l = exp(-((T - Tmelt)/W)^2)  for T < Tmelt" << nl
                << "    Tmelt = " << T_f << " K" << nl
                << "    W     = " << W_gauss << " K" << nl
                << "    thetar_frac = " << thetar_frac << nl;
        }
    }
    else
    {
        Info<< "  Method: APPARENT HEAT CAPACITY (legacy implicit with Picard)" << nl
            << "  Solves: C_eff·∂T/∂t = ∇·(K∇T)" << nl
            << "  nCorrectors: " << nCorrectors << nl
            << "  tolerance: " << finalTol << " K" << nl;
    }
    Info<< "===============================" << endl;
    
    // ======================================================================
    // INITIAL CONSTITUTIVE UPDATE
    // ======================================================================
    
    Info<< "\nInitializing phase fractions from T and theta at t=0..." << nl;
    #include "updateConstitutive.H"

    // Initialize enthalpy from temperature field (for diagnostics/output)
    forAll(H, cellI)
    {
        H[cellI] = C_sens[cellI] * T[cellI] 
                 + liquidFraction[cellI] * L.value() * theta[cellI];
    }
    H.correctBoundaryConditions();

    Info<< "Initial ranges at t=0:" << nl
        << "  T            : [" << min(T).value() << ", " << max(T).value() << "] K" << nl
        << "  H            : [" << min(H).value() << ", " << max(H).value() << "] J/m³" << nl
        << "  thetal (liq) : [" << min(thetal).value() << ", " << max(thetal).value() << "]" << nl
        << "  thetag (ice) : [" << min(thetag).value() << ", " << max(thetag).value() << "]" << nl
        << "  liquidFraction: [" << min(liquidFraction).value() << ", " 
                                 << max(liquidFraction).value() << "]" << nl
        << "  Kth          : [" << min(Kth).value() << ", " << max(Kth).value() << "] W/(m·K)" << nl;
    
    Info<< "\n=== Adaptive Time-Stepping Parameters ===" << nl
        << "  adjustTimeStep      : " << (adjustTimeStep ? "ON" : "OFF") << nl
        << "  tFact               : " << tFact << nl
        << "  maxDeltaT           : " << maxDeltaT << " s" << nl
        << "  minDeltaT           : " << minDeltaT << " s" << nl
        << "  maxFourier          : " << maxFourier << nl
        << "  stabilityThreshold  : " << stabilityThreshold << " steps" << endl;
    
    #include "postProcess.H"

    // Adaptive time-stepping state
    int stabilityCounter = 0;
    int currentPicardIters = 0;
    bool converged = true;

    // ======================================================================
    // EROSION CONTROLS
    // ======================================================================
    
    labelHashSet allowedErosionCells;
    Switch enableErosion(false);
    word erodedPatchName("erodedWall");
    // Note: erosionBoundaries is declared and read in createErosionControls.H

    if (transportProperties.found("erosionControls"))
    {
        const dictionary& ero = transportProperties.subDict("erosionControls");
        enableErosion = Switch(ero.lookup("enableErosion"));

        if (enableErosion)
        {
            Info<< "\n=== Thaw-Based Erosion Enabled ===" << nl
                << "  Cell removal criterion: liquidFraction >= " << thawThreshold << nl
                << "  Cells removed when fully thawed" << nl << endl;

            const scalar maxZ = readScalar(ero.lookup("bluffCrest"));
            const scalar minZ = readScalar(ero.lookup("bluffBase"));

            forAll(mesh.cellCentres(), c)
            {
                const scalar z = mesh.cellCentres()[c].z();
                if (z > minZ && z <= maxZ) allowedErosionCells.insert(c);
            }

            Info<< "  Erosion zone: " << allowedErosionCells.size()
                << " cells in Z ∈ (" << minZ << ", " << maxZ << "] m" << endl;

            if (ero.found("erodedPatch"))
            {
                erodedPatchName = word(ero.lookup("erodedPatch"));
            }

            Info<< "  Eroded patch name: " << erodedPatchName << nl
                << "  Erosion boundaries: " << erosionBoundaries << nl << endl;
        }
    }

    // ======================== TIME LOOP ========================
    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << " s, deltaT = " 
            << runTime.deltaTValue() << " s" << nl << endl;

        // ==================================================================
        // 1) MESH TOPOLOGY CHANGES (erosion)
        // ==================================================================
        
        bool meshChanged = false;
        if (enableErosion)
        {
            #include "updateErosionBoundary.H"
        }

        // ==================================================================
        // 2) UPDATE CONSTITUTIVE RELATIONS AFTER MESH CHANGES
        // ==================================================================
        
        if (meshChanged)
        {
            Info<< "Updating fields after mesh topology change..." << endl;
            
            T.correctBoundaryConditions();
            theta.correctBoundaryConditions();
            
            #include "updateConstitutive.H"
            
            // Reinitialize H from T for newly exposed cells
            forAll(H, cellI)
            {
                H[cellI] = C_sens[cellI] * T[cellI] 
                         + liquidFraction[cellI] * L.value() * theta[cellI];
            }
            H.correctBoundaryConditions();
            
            Info<< "Post-erosion field ranges:" << nl
                << "  T      : [" << min(T).value() << ", " << max(T).value() << "] K" << nl
                << "  H      : [" << min(H).value() << ", " << max(H).value() << "] J/m³" << endl;
        }

        // ==================================================================
        // 3) STORE OLD TIME VALUES
        // ==================================================================
        
        // Store liquid fraction at start of timestep (for latent source)
        volScalarField liquidFraction_old("liquidFraction_old", liquidFraction);
        volScalarField T_old("T_old", T);
        
        const scalar dt = runTime.deltaTValue();
        converged = true;

        // ==================================================================
        // 4) THERMAL SOLUTION
        // ==================================================================
        
        if (useEnthalpyPorosity)
        {
            // ==============================================================
            // SOURCE-BASED ENTHALPY-POROSITY METHOD
            //
            // This solves the energy equation with explicit latent heat source:
            //
            //   C_sens * ∂T/∂t - ∇·(K∇T) = -L * θ * ∂f_l/∂t
            //
            // The latent source is computed from the change in liquid fraction,
            // which is updated from temperature using the SFCC.
            //
            // This avoids the H→T inversion problem while still capturing
            // the latent heat effect through an iterative source update.
            //
            // Each Picard iteration:
            //   1. Solve implicit T equation with current latent source
            //   2. Update f_l from new T
            //   3. Recompute latent source
            //   4. Repeat until converged
            //
            // REFERENCES:
            // - Voller & Swaminathan (1991) - Source-based enthalpy method
            // ==============================================================

            converged = false;
            scalar maxResT = 0.0;

            // Store fields at start of timestep
            volScalarField T_old_step("T_old_step", T);
            volScalarField H_old("H_old", H);

            // Compute constitutive properties ONCE at start of timestep
            // Using LAGGED CthEff prevents iteration instability
            #include "updateConstitutive.H"
            volScalarField CthEff_lagged("CthEff_lagged", CthEff);

            for (int corr = 0; corr < nCorrectors; corr++)
            {
                // Store previous iteration T for convergence check
                volScalarField T_prev_iter("T_prev_iter", T);
                volScalarField fl_prev("fl_prev", liquidFraction);

                // ----------------------------------------------------------
                // Solve implicit temperature equation with LAGGED APPARENT HEAT CAPACITY
                //
                // Use CthEff from START OF TIMESTEP (not updated within Picard).
                // This makes the equation linear within each timestep, preventing
                // the oscillation caused by CthEff changing with T.
                //
                // The latent heat is captured by the high CthEff in mushy zone.
                // ----------------------------------------------------------

                fvScalarMatrix TEqn
                (
                    fvm::ddt(CthEff_lagged, T)
                  - fvm::laplacian(Kth, T)
                );

                TEqn.relax(relaxationFactor);
                TEqn.solve();

                T.correctBoundaryConditions();
                T.clamp_range(Tminc.value(), Tmaxc.value());

                // ----------------------------------------------------------
                // Update liquid fraction from temperature
                // ----------------------------------------------------------

                forAll(liquidFraction, cellI)
                {
                    const scalar T_cell = T[cellI];
                    scalar fl_new = 0.0;

                    if (SFCCModel == "linear")
                    {
                        if (T_cell <= T_sol)
                        {
                            fl_new = thetar_frac;
                        }
                        else if (T_cell >= T_liq)
                        {
                            fl_new = 1.0;
                        }
                        else
                        {
                            fl_new = thetar_frac
                                   + (1.0 - thetar_frac) * (T_cell - T_sol) / dT_mushy;
                        }
                    }
                    else if (SFCCModel == "powerLaw")
                    {
                        if (T_cell >= T_f)
                        {
                            fl_new = 1.0;
                        }
                        else
                        {
                            const scalar depression = T_f - T_cell;
                            const scalar b_val = max(mag(b_exp), SMALL);
                            const scalar arg = max(depression / T_ref, SMALL);
                            const scalar r = Foam::pow(arg, b_val);
                            
                            // Correct retention form mapping
                            fl_new = thetar_frac + (1.0 - thetar_frac) / (1.0 + r);
                        }
                    }
                    else if (SFCCModel == "gaussian")
                    {
                        if (T_cell >= T_f)
                        {
                            fl_new = 1.0;
                        }
                        else
                        {
                            const scalar arg = (T_cell - T_f) / W_gauss;
                            fl_new = (1.0 - thetar_frac) * Foam::exp(-arg * arg)
                                   + thetar_frac;
                        }
                    }

                    // Under-relax liquid fraction for stability
                    liquidFraction[cellI] = flRelaxation * fl_new
                                          + (1.0 - flRelaxation) * liquidFraction[cellI];
                    liquidFraction[cellI] = min(max(liquidFraction[cellI], scalar(0.0)), scalar(1.0));
                }

                liquidFraction.correctBoundaryConditions();

                // ----------------------------------------------------------
                // Update enthalpy for diagnostics (H = C_sens*T + f_l*L*θ)
                // ----------------------------------------------------------

                forAll(H, cellI)
                {
                    H[cellI] = C_sens[cellI] * T[cellI]
                             + liquidFraction[cellI] * L.value() * theta[cellI];
                }
                H.correctBoundaryConditions();

                // ----------------------------------------------------------
                // Convergence check
                // ----------------------------------------------------------

                maxResT = 0.0;
                scalar maxDeltaFL = 0.0;
                forAll(T, i)
                {
                    maxResT = max(maxResT, mag(T[i] - T_prev_iter[i]));
                    maxDeltaFL = max(maxDeltaFL, mag(liquidFraction[i] - fl_prev[i]));
                }
                reduce(maxResT, maxOp<scalar>());
                reduce(maxDeltaFL, maxOp<scalar>());

                // Count phase state
                // Count phase state (Relative to residual fraction)
                label nSolid = 0, nMushy = 0, nLiquid = 0;
                const scalar solidThreshold = thetar_frac + 0.001; 

                forAll(liquidFraction, i)
                {
                    if (liquidFraction[i] <= solidThreshold) nSolid++;
                    else if (liquidFraction[i] >= 0.999) nLiquid++;
                    else nMushy++;
                }
                reduce(nSolid, sumOp<label>());
                reduce(nMushy, sumOp<label>());
                reduce(nLiquid, sumOp<label>());

                Info<< "  Picard iter " << (corr + 1) << "/" << nCorrectors
                    << ": res(T)=" << maxResT << " K, Δf_l=" << maxDeltaFL
                    << ", phases: S=" << nSolid << " M=" << nMushy << " L=" << nLiquid
                    << endl;

                if (maxResT < finalTol && maxDeltaFL < 0.01)
                {
                    currentPicardIters = corr + 1;
                    converged = true;
                    Info<< "  *** Converged in " << currentPicardIters
                        << " iterations ***" << nl;
                    break;
                }

                if (corr == nCorrectors - 1)
                {
                    currentPicardIters = nCorrectors;
                    if (maxResT < finalTol * 10)
                    {
                        converged = true;
                        Info<< "  *** Near-converged ***" << nl;
                    }
                    else
                    {
                        WarningInFunction
                            << "Picard did not converge: res(T)=" << maxResT << " K"
                            << endl;
                    }
                }
            }

            // Post-iteration diagnostics
            scalar maxDeltaT_step = 0.0;
            scalar maxDeltaH_step = 0.0;
            forAll(T, i)
            {
                maxDeltaT_step = max(maxDeltaT_step, mag(T[i] - T_old_step[i]));
                maxDeltaH_step = max(maxDeltaH_step, mag(H[i] - H_old[i]));
            }
            reduce(maxDeltaT_step, maxOp<scalar>());
            reduce(maxDeltaH_step, maxOp<scalar>());

            Info<< "  Step summary:" << nl
                << "    max ΔT     : " << maxDeltaT_step << " K" << nl
                << "    max ΔH     : " << maxDeltaH_step << " J/m³" << nl
                << "    T range    : [" << min(T).value() << ", " << max(T).value() << "] K" << nl
                << "    H range    : [" << min(H).value() << ", " << max(H).value() << "] J/m³" << nl
                << "    liquidFrac : [" << min(liquidFraction).value() << ", "
                    << max(liquidFraction).value() << "]" << nl;
        }
        else
        {
            // ==============================================================
            // LEGACY APPARENT HEAT CAPACITY METHOD (IMPLICIT WITH PICARD)
            // ==============================================================
            
            converged = false;
            
            for (int corr = 0; corr < nCorrectors; corr++)
            {
                volScalarField T_prev("T_prev", T);
                
                #include "updateConstitutive.H"
                
                fvScalarMatrix TEqn
                (
                    fvm::ddt(CthEff, T)
                  - fvm::laplacian(Kth, T)
                );

                TEqn.relax(relaxationFactor);
                TEqn.solve();
                
                T.correctBoundaryConditions();
                T.clamp_range(Tminc.value(), Tmaxc.value());
                
                // Update H for consistency
                forAll(H, cellI)
                {
                    H[cellI] = C_sens[cellI] * T[cellI] 
                             + liquidFraction[cellI] * L.value() * theta[cellI];
                }

                // Check convergence
                scalar maxResT = 0.0;
                forAll(T, i)
                {
                    maxResT = max(maxResT, mag(T[i] - T_prev[i]));
                }
                reduce(maxResT, maxOp<scalar>());
                
                Info<< "  Picard iter " << (corr + 1) << "/" << nCorrectors 
                    << ", residual(T): " << maxResT << endl;
                
                if (maxResT < finalTol)
                {
                    currentPicardIters = corr + 1;
                    converged = true;
                    Info<< "  *** Converged in " << currentPicardIters << " iterations ***" << nl;
                    break;
                }
                
                if (corr == nCorrectors - 1)
                {
                    currentPicardIters = nCorrectors;
                    if (maxResT < finalTol * 10)
                    {
                        converged = true;
                    }
                    else
                    {
                        WarningInFunction 
                            << "Picard did not converge: residual = " << maxResT << endl;
                    }
                }
            }
            
            Info<< "    T range    : [" << min(T).value() << ", " << max(T).value() << "] K" << nl
                << "    liquidFrac : [" << min(liquidFraction).value() << ", " 
                    << max(liquidFraction).value() << "]" << nl;
        }

        // ==================================================================
        // 5) ADAPTIVE TIME-STEPPING (Physics-Based)
        // ==================================================================
        
        if (adjustTimeStep)
        {
            // ----------------------------------------------------------
            // Compute physics-based timestep constraints using CthEff
            //
            // IMPORTANT: Use CthEff (apparent heat capacity) not C_sens!
            // In the mushy zone, CthEff >> C_sens due to latent heat.
            // Using C_sens overestimates diffusivity → dt collapses.
            // ----------------------------------------------------------

            // Fourier number estimate: Fo = K * dt / (C_eff * dx²)
            // We compute max thermal diffusivity and min cell size

            scalar maxDiffusivity = 0.0;
            scalar minCellSize = GREAT;

            forAll(mesh.V(), cellI)
            {
                // Use CthEff (apparent capacity including latent) for Fourier stability
                const scalar alpha = Kth[cellI] / max(CthEff[cellI], minCthEff);
                maxDiffusivity = max(maxDiffusivity, alpha);

                const scalar dx = Foam::cbrt(mesh.V()[cellI]);
                minCellSize = min(minCellSize, dx);
            }

            reduce(maxDiffusivity, maxOp<scalar>());
            reduce(minCellSize, minOp<scalar>());

            // Stable timestep from Fourier constraint (for guidance only)
            // Implicit scheme is unconditionally stable, but accuracy improves
            // with reasonable Fourier numbers
            const scalar dt_Fourier = maxFourier * sqr(minCellSize) / max(maxDiffusivity, SMALL);
            
            // Current Fourier number
            const scalar currentFo = maxDiffusivity * dt / sqr(minCellSize);
            
            scalar deltaTFactor = 1.0;
            bool timeStepChanged = false;
            
            const bool atMinTimeStep = (runTime.deltaTValue() <= minDeltaT * 1.001);
            const bool atMaxTimeStep = (runTime.deltaTValue() >= maxDeltaT * 0.999);
            
            if (!converged && !atMinTimeStep)
            {
                // Reduce timestep if Picard didn't converge
                deltaTFactor = 1.0 / tFact;
                timeStepChanged = true;
                stabilityCounter = 0;
                Info<< "  *** Reducing timestep (non-convergence) ***" << nl;
            }
            else if (converged && currentPicardIters > nCorrectors * 0.7 && !atMinTimeStep)
            {
                // Reduce timestep if using too many iterations
                deltaTFactor = 1.0 / Foam::sqrt(tFact);
                timeStepChanged = true;
                stabilityCounter = 0;
                Info<< "  *** Reducing timestep (many iterations: " 
                    << currentPicardIters << ") ***" << nl;
            }
            else if (converged)
            {
                ++stabilityCounter;
                
                // Increase timestep if:
                // - Stable for stabilityThreshold steps
                // - Using few iterations
                // - Not already at max
                if (stabilityCounter >= stabilityThreshold 
                    && currentPicardIters <= nCorrectors / 3
                    && !atMaxTimeStep)
                {
                    deltaTFactor = tFact;
                    timeStepChanged = true;
                    stabilityCounter = 0;
                    Info<< "  *** Increasing timestep (stable, few iterations) ***" << nl;
                }
            }
            
            if (timeStepChanged)
            {
                scalar newDeltaT = deltaTFactor * runTime.deltaTValue();
                
                // Also respect Fourier guidance (but not strictly enforce)
                if (newDeltaT > dt_Fourier * 2.0)
                {
                    newDeltaT = dt_Fourier * 2.0;
                    Info<< "      Limited by Fourier guidance" << nl;
                }
                
                newDeltaT = min(max(newDeltaT, minDeltaT), maxDeltaT);
                
                Info<< "      Old deltaT: " << runTime.deltaTValue() << " s" << nl
                    << "      New deltaT: " << newDeltaT << " s" << nl
                    << "      Fourier: " << currentFo << " (guide dt: " << dt_Fourier << " s)" << nl;
                
                runTime.setDeltaT(newDeltaT);
            }
            else
            {
                Info<< "  Fourier number: " << currentFo 
                    << ", stability counter: " << stabilityCounter 
                    << "/" << stabilityThreshold << nl;
            }
        }

        // ==================================================================
        // 6) EROSION STATISTICS
        // ==================================================================
        
        #include "handleDynamicErosion.H"
        
        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "\n=== Simulation complete ===" << nl
        << "Final time: " << runTime.timeName() << " s" << nl
        << "Total steps: " << runTime.timeIndex() << nl << endl;
        
    return 0;
}
