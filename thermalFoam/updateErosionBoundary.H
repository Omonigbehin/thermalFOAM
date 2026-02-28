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
    updateErosionBoundary.H - Simple Thaw-Based Erosion Model (PARALLEL-SAFE)

    Simplified erosion model: cells are removed when fully thawed.

    This replaces the D*-based model from Omonigbehin et al. (2025) with a
    simpler phase-change-driven approach:
      - Cells are removed when liquidFraction >= thawThreshold
      - No erosion rate law - immediate removal upon reaching threshold
      - Niche morphology tracking (depth and volume) is retained

    REQUIRED INPUTS:
      From createErosionControls.H:
        thawThreshold, bluffCrest, bluffBase,
        nicheTrackingZ, nicheTrackingTol

      From thermalFoam.C:
        enableErosion, erosionBoundaries, erodedPatchName, allowedErosionCells

      From createFields.H:
        T, liquidFraction, erosionMask
\*---------------------------------------------------------------------------*/

#include "removeCells.H"
#include "polyTopoChange.H"
#include "syncTools.H"

Info<< "--- Thaw-Based Erosion Check ---" << nl;

// ==========================================================================
// STATIC VARIABLES FOR CUMULATIVE TRACKING
// ==========================================================================

static scalar cumulativeErodedVolume = 0.0;
static label cumulativeErodedCells = 0;
static scalar initialNicheX = -GREAT;  // Track initial bluff position at niche tracking Z
static bool erosionFileInitialized = false;

// ==========================================================================
// INITIALIZE BLUFF POSITION TRACKING AT NICHE Z (first call only)
// ALL PROCESSORS MUST EXECUTE THIS BLOCK TOGETHER
// ==========================================================================

if (initialNicheX < -1e10)
{
    // Find initial bluff face position at the niche tracking elevation
    // We search for the MAXIMUM X at the tracking Z (assuming +X points inland)
    // LOCAL maximum first, then global reduce

    scalar localNicheX = -GREAT;

    forAll(erosionBoundaries, i)
    {
        const word& patchName = erosionBoundaries[i];
        const label pID = mesh.boundaryMesh().findPatchID(patchName);
        if (pID == -1) continue;

        const polyPatch& pp = mesh.boundaryMesh()[pID];

        forAll(pp.faceCentres(), fI)
        {
            const vector& fc = pp.faceCentres()[fI];

            // Only consider faces near the niche tracking elevation
            if (mag(fc.z() - nicheTrackingZ) <= nicheTrackingTol)
            {
                localNicheX = max(localNicheX, fc.x());
            }
        }
    }

    // GLOBAL reduction - ALL processors must call this
    initialNicheX = localNicheX;
    reduce(initialNicheX, maxOp<scalar>());

    if (initialNicheX > -1e10)
    {
        if (Pstream::master())
        {
            Info<< "  Initial bluff X at Z=" << nicheTrackingZ << " m: "
                << initialNicheX << " m" << nl;
        }
    }
    else
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Could not find initial bluff position at Z = " << nicheTrackingZ
                << " ± " << nicheTrackingTol << " m" << nl
                << "Check that erosionBoundaries patches have faces at this elevation." << nl
                << "Niche depth tracking may be inaccurate." << endl;
        }

        // Fallback: use any face in erosion zone
        scalar localFallbackX = -GREAT;

        forAll(erosionBoundaries, i)
        {
            const word& patchName = erosionBoundaries[i];
            const label pID = mesh.boundaryMesh().findPatchID(patchName);
            if (pID == -1) continue;

            const polyPatch& pp = mesh.boundaryMesh()[pID];

            forAll(pp.faceCentres(), fI)
            {
                const vector& fc = pp.faceCentres()[fI];
                if (fc.z() > bluffBase && fc.z() <= bluffCrest)
                {
                    localFallbackX = max(localFallbackX, fc.x());
                }
            }
        }

        // GLOBAL reduction for fallback
        initialNicheX = localFallbackX;
        reduce(initialNicheX, maxOp<scalar>());

        if (Pstream::master())
        {
            Info<< "  Fallback: Initial bluff X (any Z in erosion zone): "
                << initialNicheX << " m" << nl;
        }
    }
}

// ==========================================================================
// INITIALIZE EROSION LOG FILE (first call only, master only)
// ==========================================================================

if (!erosionFileInitialized)
{
    if (Pstream::master())
    {
        mkDir("postProcessing");
        OFstream erosionHeader("postProcessing/erosionHistory.dat");
        erosionHeader
            << "# Thaw-based erosion history" << nl
            << "# Cells removed when liquidFraction >= " << thawThreshold << nl
            << "# Niche depth tracked at Z = " << nicheTrackingZ
            << " ± " << nicheTrackingTol << " m" << nl
            << "# Time[s]  Volume[m³]  Cells[-]  NicheDepth[m]  "
            << "StepVolume[m³]  StepCells[-]  MeanLiquidFrac[-]  MaxLiquidFrac[-]" << nl;
    }
    // ALL processors set the flag
    erosionFileInitialized = true;
}

// ==========================================================================
// COLLECT EROSION CANDIDATE CELLS
// ==========================================================================

labelHashSet candidateCells;

// Source 1: Original erosion boundaries
forAll(erosionBoundaries, i)
{
    const word& patchName = erosionBoundaries[i];
    const label pID = mesh.boundaryMesh().findPatchID(patchName);
    if (pID == -1) continue;

    const polyPatch& pp = mesh.boundaryMesh()[pID];
    const labelList& faceCells = pp.faceCells();

    forAll(faceCells, fI)
    {
        const label c = faceCells[fI];
        if (allowedErosionCells.found(c))
        {
            candidateCells.insert(c);
        }
    }
}

// Source 2: Eroded wall patch (progressive inland erosion)
const label erodedPatchID = mesh.boundaryMesh().findPatchID(erodedPatchName);
if (erodedPatchID != -1)
{
    const polyPatch& pp = mesh.boundaryMesh()[erodedPatchID];
    const labelList& faceCells = pp.faceCells();

    forAll(faceCells, fI)
    {
        const label c = faceCells[fI];
        if (allowedErosionCells.found(c))
        {
            candidateCells.insert(c);
        }
    }

    if (pp.size() > 0)
    {
        Info<< "  Eroded wall: " << pp.size() << " faces" << nl;
    }
}

// Report candidate count (local, for debugging)
label localCandidates = candidateCells.size();
label globalCandidates = localCandidates;
reduce(globalCandidates, sumOp<label>());
Info<< "  Candidate cells: " << globalCandidates << " (local: " << localCandidates << ")" << nl;

// ==========================================================================
// SIMPLE THAW-BASED REMOVAL CRITERION
// ==========================================================================

DynamicList<label> cellsToRemove;
scalar maxLiquidFrac = 0.0;
scalar sumLiquidFrac = 0.0;
label nThawedCells = 0;

forAllConstIter(labelHashSet, candidateCells, iter)
{
    const label c = iter.key();
    const scalar lf_cell = liquidFraction[c];

    // Track statistics
    maxLiquidFrac = max(maxLiquidFrac, lf_cell);
    sumLiquidFrac += lf_cell;

    // ======================================================================
    // SIMPLE CRITERION: Remove cell if thawed beyond threshold
    // ======================================================================

    if (lf_cell >= thawThreshold)
    {
        cellsToRemove.append(c);
        erosionMask[c] = 1.0;
        nThawedCells++;
    }
}

// ==========================================================================
// GLOBAL STATISTICS - ALL processors participate
// ==========================================================================

scalar globalMaxLiquidFrac = maxLiquidFrac;
scalar globalSumLiquidFrac = sumLiquidFrac;
label globalThawedCells = nThawedCells;

reduce(globalMaxLiquidFrac, maxOp<scalar>());
reduce(globalSumLiquidFrac, sumOp<scalar>());
reduce(globalThawedCells, sumOp<label>());

const scalar meanLiquidFrac = (globalCandidates > 0)
                            ? globalSumLiquidFrac / globalCandidates
                            : 0.0;

Info<< "  LiquidFraction statistics: max=" << globalMaxLiquidFrac
    << ", mean=" << meanLiquidFrac
    << " (" << globalThawedCells << " cells >= threshold " << thawThreshold << ")" << nl;

// ==========================================================================
// SYNC CELLS TO REMOVE ACROSS PROCESSORS
// ==========================================================================

label localRemoveCount = cellsToRemove.size();
label globalRemoveCount = localRemoveCount;
reduce(globalRemoveCount, sumOp<label>());

Info<< "  Cells to remove (thawed): " << globalRemoveCount
    << " (local: " << localRemoveCount << ")" << nl;

// ==========================================================================
// PERFORM MESH MODIFICATION IF NEEDED
// ALL PROCESSORS MUST EXECUTE THIS BLOCK TOGETHER (even if they have 0 cells)
// ==========================================================================

if (globalRemoveCount > 0)
{
    // Volume tracking (LOCAL first, then reduce)
    scalar localStepVolume = 0.0;
    forAll(cellsToRemove, i)
    {
        localStepVolume += mesh.V()[cellsToRemove[i]];
    }

    // GLOBAL reduction for volume
    scalar stepVolume = localStepVolume;
    reduce(stepVolume, sumOp<scalar>());

    cumulativeErodedVolume += stepVolume;
    cumulativeErodedCells += globalRemoveCount;

    // ======================================================================
    // MESH TOPOLOGY CHANGE - COLLECTIVE OPERATION
    // ======================================================================

    removeCells cutter(mesh);

    // Find/create exposed patch
    label exposedPatchID = mesh.boundaryMesh().findPatchID(erodedPatchName);
    if (exposedPatchID == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch '" << erodedPatchName << "' for exposed faces"
            << abort(FatalError);
    }

    // Get cells and exposed faces for THIS processor
    const labelList cellsLab(cellsToRemove);
    const labelList facesToExpose = cutter.getExposedFaces(cellsLab);
    labelList patchIDs(facesToExpose.size(), exposedPatchID);

    // Perform mesh modification - ALL processors must call this
    polyTopoChange meshMod(mesh);
    cutter.setRefinement(cellsLab, facesToExpose, patchIDs, meshMod);

    // changeMesh is a COLLECTIVE operation - ALL processors must call
    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    cutter.updateMesh(morphMap());
    mesh.updateMesh(morphMap());
    mesh.update();

    // ======================================================================
    // REBUILD allowedErosionCells (indices changed)
    // ======================================================================

    allowedErosionCells.clear();
    forAll(mesh.C(), cellI)
    {
        const scalar zCoord = mesh.C()[cellI].z();
        if (zCoord > bluffBase && zCoord <= bluffCrest)
        {
            allowedErosionCells.insert(cellI);
        }
    }

    // ======================================================================
    // FIX FIELDS AFTER TOPOLOGY CHANGE
    // ======================================================================

    forAll(Kth, cellI)
    {
        if (!std::isfinite(Kth[cellI]) || Kth[cellI] < 1e-6)
        {
            Kth[cellI] = 1.5;  // Default thermal conductivity
        }
    }
    Kth.correctBoundaryConditions();

    forAll(T, cellI)
    {
        if (!std::isfinite(T[cellI]))
        {
            T[cellI] = Tmelt.value();
            WarningInFunction
                << "NaN/Inf in T at cell " << cellI << ", reset to Tmelt" << endl;
        }
        T[cellI] = min(max(T[cellI], Tminc.value()), Tmaxc.value());
    }

    T.correctBoundaryConditions();
    theta.correctBoundaryConditions();

    // ======================================================================
    // COMPUTE NICHE DEPTH AT SPECIFIED Z-COORDINATE
    // ======================================================================

    scalar currentNicheX = -GREAT;
    label nFacesAtTrackingZ = 0;

    const label finalErodedPatchID = mesh.boundaryMesh().findPatchID(erodedPatchName);
    if (finalErodedPatchID != -1)
    {
        const polyPatch& erodedPatch = mesh.boundaryMesh()[finalErodedPatchID];

        // Search only faces near the niche tracking elevation
        forAll(erodedPatch.faceCentres(), fI)
        {
            const vector& fc = erodedPatch.faceCentres()[fI];

            if (mag(fc.z() - nicheTrackingZ) <= nicheTrackingTol)
            {
                currentNicheX = max(currentNicheX, fc.x());
                nFacesAtTrackingZ++;
            }
        }
    }

    // Global reductions for niche tracking
    reduce(currentNicheX, maxOp<scalar>());
    reduce(nFacesAtTrackingZ, sumOp<label>());

    // Niche depth = inland penetration from initial position
    // Positive value indicates erosion into the bluff
    scalar nicheDepth = 0.0;

    if (initialNicheX > -1e10 && currentNicheX > -1e10)
    {
        nicheDepth = currentNicheX - initialNicheX;

        // If negative, the geometry convention might be reversed
        // (i.e., +X points seaward instead of inland)
        if (nicheDepth < 0)
        {
            nicheDepth = -nicheDepth;  // Use absolute value
        }
    }

    // Warn if no faces found at tracking elevation (master only)
    if (nFacesAtTrackingZ == 0 && finalErodedPatchID != -1)
    {
        static bool warnedOnce = false;
        if (!warnedOnce && Pstream::master())
        {
            Info<< "  NOTE: No eroded faces at Z=" << nicheTrackingZ
                << " m yet. Niche depth = 0 until erosion reaches this level." << nl;
            warnedOnce = true;
        }
    }

    // ======================================================================
    // VOLUME REPORTING (master only)
    // ======================================================================

    const scalar postErosionVol = gSum(mesh.V());

    if (Pstream::master())
    {
        Info<< "  Step eroded volume : " << stepVolume << " m³" << nl
            << "  Cumulative volume  : " << cumulativeErodedVolume << " m³" << nl
            << "  Cumulative cells   : " << cumulativeErodedCells << nl
            << "  Niche depth (Z=" << nicheTrackingZ << "m): "
            << nicheDepth*1000 << " mm" << nl
            << "  Faces at tracking Z: " << nFacesAtTrackingZ << nl
            << "  Remaining volume   : " << postErosionVol << " m³" << endl;
    }

    // ======================================================================
    // WRITE TO CUMULATIVE LOG FILE (master only)
    // ======================================================================

    if (Pstream::master())
    {
        OFstream erosionLog
        (
            "postProcessing/erosionHistory.dat",
            IOstreamOption(),
            IOstreamOption::APPEND
        );

        erosionLog
            << runTime.value() << "  "
            << cumulativeErodedVolume << "  "
            << cumulativeErodedCells << "  "
            << nicheDepth << "  "
            << stepVolume << "  "
            << globalRemoveCount << "  "
            << meanLiquidFrac << "  "
            << globalMaxLiquidFrac << nl;
    }

    meshChanged = true;
}
else
{
    Info<< "  No cells removed this timestep" << nl;
}

Info<< "--- End thaw-based erosion ---" << nl;
