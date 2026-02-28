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
    handleDynamicErosion.H - Erosion Statistics and Safety Monitoring
    (PARALLEL-SAFE VERSION)

    This file handles:
    1. Runtime safety checks (excessive erosion warning)
    2. Output-time diagnostics

    Note: Cumulative erosion logging is now handled in updateErosionBoundary.H
    to ensure logging occurs exactly when cells are removed.

    PARALLEL FIX: All reduce() operations are now UNCONDITIONAL to prevent
    deadlocks when runTime.outputTime() differs across processors.
\*---------------------------------------------------------------------------*/

if (enableErosion)
{
    // ======================================================================
    // TRACK INITIAL DOMAIN VOLUME (first call only)
    // ======================================================================

    static scalar initialDomainVolume = -1.0;
    if (initialDomainVolume < 0)
    {
        initialDomainVolume = gSum(mesh.V());
        if (Pstream::master())
        {
            Info<< "Initial domain volume: " << initialDomainVolume << " m³" << endl;
        }
    }

    // ======================================================================
    // SAFETY CHECK: Warn if excessive erosion
    // ======================================================================

    const scalar erosionWarningThreshold =
        transportProperties.lookupOrDefault<scalar>("erosionWarningThreshold", 0.5);

    const scalar currentVol = gSum(mesh.V());
    const scalar volumeFraction = currentVol / initialDomainVolume;

    if (volumeFraction < erosionWarningThreshold)
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Excessive erosion detected!" << nl
                << "  Current volume: " << currentVol << " m³" << nl
                << "  Initial volume: " << initialDomainVolume << " m³" << nl
                << "  Remaining: " << (volumeFraction * 100.0) << "%" << nl
                << "  Threshold: " << (erosionWarningThreshold * 100.0) << "%" << endl;
        }
    }

    // ======================================================================
    // OUTPUT-TIME DIAGNOSTICS
    // PARALLEL FIX: Compute reductions UNCONDITIONALLY, only gate the output
    // ======================================================================

    // Count currently marked erosion cells (ALL processors compute this)
    label nMarkedCells = 0;
    forAll(erosionMask, cellI)
    {
        if (erosionMask[cellI] > 0.5)
        {
            nMarkedCells++;
        }
    }

    // Liquid fraction statistics for candidate cells (ALL processors compute this)
    scalar maxLiquidFrac = 0.0;
    scalar sumLiquidFrac = 0.0;
    label nThawingCells = 0;

    forAll(liquidFraction, cellI)
    {
        if (liquidFraction[cellI] > SMALL)
        {
            maxLiquidFrac = max(maxLiquidFrac, liquidFraction[cellI]);
            sumLiquidFrac += liquidFraction[cellI];
            nThawingCells++;
        }
    }

    // CRITICAL: ALL processors must call reduce() together - NO conditional!
    reduce(nMarkedCells, sumOp<label>());
    reduce(maxLiquidFrac, maxOp<scalar>());
    reduce(sumLiquidFrac, sumOp<scalar>());
    reduce(nThawingCells, sumOp<label>());

    const scalar meanLiquidFrac = (nThawingCells > 0)
                                ? sumLiquidFrac / nThawingCells
                                : 0.0;

    // Only PRINT at output times, but reductions already done
    if (runTime.outputTime() && Pstream::master())
    {
        Info<< "=== Erosion Status at t = " << runTime.timeName() << " ===" << nl
            << "  Domain volume    : " << currentVol << " m³ ("
            << (volumeFraction * 100.0) << "% remaining)" << nl
            << "  Erosion mask     : " << nMarkedCells << " cells marked" << nl
            << "  Liquid fraction  : max = " << maxLiquidFrac
            << ", mean = " << meanLiquidFrac << nl
            << "  Thawing cells    : " << nThawingCells << nl
            << "  Thaw threshold   : " << thawThreshold << nl
            << "=======================================" << endl;
    }
}
