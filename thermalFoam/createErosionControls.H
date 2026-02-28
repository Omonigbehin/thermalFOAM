/*---------------------------------------------------------------------------*\
  Copyright (C) 2025 Olorunfemi Omonigbehin, INRS-ETE

    createErosionControls.H - Simple Thaw-Based Erosion Parameters
      - Cells are removed when liquidFraction >= thawThreshold
      - No erosion rate law - immediate removal upon reaching threshold
      - Niche morphology tracking (depth and volume) is retained

    This is an helper file that should be included in createFields.H 
\*---------------------------------------------------------------------------*/

// ==========================================================================
// THAW-BASED EROSION PARAMETERS
// ==========================================================================

// Erosion boundary patch names (used for originalPatchID tracking)
wordList erosionBoundaries;

// Thaw threshold for cell removal (0.0 to 1.0)
// 1.0 = fully thawed (all ice melted)
// 0.5 = half melted
// 0.0 = remove as soon as any melting begins
scalar thawThreshold = 1.0; //default value

// Geometry bounds for erosion zone (Z-coordinates)
scalar bluffCrest = 0.20;
scalar bluffBase = 0.11;

// ==========================================================================
// NICHE DEPTH TRACKING PARAMETERS
// ==========================================================================

scalar nicheTrackingZ = 0.135;      // Z-coordinate to track niche depth [m]
scalar nicheTrackingTol = 0.005;   // Tolerance band around Z [m]

// Read from dictionary if erosionControls exists
if (transportProperties.found("erosionControls"))
{
    const dictionary& erosionControlsDict = transportProperties.subDict("erosionControls");

    Info<< "\nReading thaw-based erosion parameters..." << endl;

    // ======================================================================
    // THAW THRESHOLD
    // ======================================================================

    thawThreshold = erosionControlsDict.lookupOrDefault<scalar>("thawThreshold", 1.0);

    // Ensure valid range [0, 1]
    thawThreshold = max(scalar(0.0), min(thawThreshold, scalar(1.0)));

    Info<< "  Thaw threshold: liquidFraction >= " << thawThreshold << nl;

    if (thawThreshold < 1.0)
    {
        Info<< "  NOTE: Cells will be removed BEFORE fully thawed" << nl;
    }
    else
    {
        Info<< "  Cells removed only when FULLY thawed (liquidFraction = 1.0)" << nl;
    }

    // ======================================================================
    // GEOMETRY BOUNDS
    // ======================================================================

    bluffCrest = erosionControlsDict.lookupOrDefault<scalar>("bluffCrest", 0.20);
    bluffBase = erosionControlsDict.lookupOrDefault<scalar>("bluffBase", 0.11);

    Info<< "  Erosion zone: Z ∈ (" << bluffBase << ", " << bluffCrest << "] m" << endl;

    // ======================================================================
    // NICHE DEPTH TRACKING PARAMETERS
    // ======================================================================

    // Default: track at 90% of the way from bluffBase to bluffCrest
    // This is typically where maximum niche penetration occurs
    const scalar defaultNicheZ = 0.9 * bluffCrest + 0.1 * bluffBase;

    nicheTrackingZ = erosionControlsDict.lookupOrDefault<scalar>("nicheTrackingZ", defaultNicheZ);
    nicheTrackingTol = erosionControlsDict.lookupOrDefault<scalar>("nicheTrackingTol", 0.005);

    Info<< "  Niche depth tracking:" << nl
        << "    Z-coordinate      = " << nicheTrackingZ << " m" << nl
        << "    Tolerance band    = ±" << nicheTrackingTol*1000 << " mm" << endl;

    // Validate niche tracking Z is within erosion zone
    if (nicheTrackingZ < bluffBase || nicheTrackingZ > bluffCrest)
    {
        WarningInFunction
            << "nicheTrackingZ = " << nicheTrackingZ << " m is outside erosion zone" << nl
            << "Erosion zone: Z ∈ (" << bluffBase << ", " << bluffCrest << "] m" << nl
            << "Niche depth tracking may not capture erosion progress." << endl;
    }

    // ======================================================================
    // EROSION BOUNDARY PATCHES
    // ======================================================================

    if (erosionControlsDict.found("erosionBoundaries"))
    {
        erosionBoundaries = erosionControlsDict.lookup("erosionBoundaries");
        Info<< "  Erosion boundary patches: " << erosionBoundaries << endl;
    }

    // ======================================================================
    // these are palceholders for future extensions
    // ======================================================================

    if (erosionControlsDict.found("waveHeight") ||
        erosionControlsDict.found("wavePeriod") ||
        erosionControlsDict.found("nearshoreSlope") ||
        erosionControlsDict.found("waterDepthAtBluff"))
    {
        WarningInFunction
            << "D* wave parameters found but are NO LONGER USED." << nl
            << "The thaw-based erosion model removes cells when liquidFraction >= thawThreshold." << nl
            << "You can safely remove waveHeight, wavePeriod, nearshoreSlope, waterDepthAtBluff." << endl;
    }

    if (erosionControlsDict.found("Cp_permafrost") ||
        erosionControlsDict.found("Lp_permafrost") ||
        erosionControlsDict.found("waterTemperature"))
    {
        WarningInFunction
            << "D* thermal parameters found but are NO LONGER USED." << nl
            << "The thaw-based erosion model does not use Stefan number calculations." << nl
            << "You can safely remove Cp_permafrost, Lp_permafrost, waterTemperature." << endl;
    }

    if (erosionControlsDict.found("useLiquidFractionGate") ||
        erosionControlsDict.found("lfOnset") ||
        erosionControlsDict.found("lfFull") ||
        erosionControlsDict.found("erosionTempOffset") ||
        erosionControlsDict.found("erosionGateWidth"))
    {
        WarningInFunction
            << "Smoothstep gate parameters found but are NO LONGER USED." << nl
            << "The thaw-based erosion model uses a simple threshold: thawThreshold." << nl
            << "Set thawThreshold to control when cells are removed (default 1.0 = fully thawed)." << endl;
    }
}

Info<< "Thaw-based erosion parameters initialized.\n" << endl;

