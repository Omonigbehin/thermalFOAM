/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2025 Femi Omonigbehin, INRS-ETE (wave extensions)
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "waveExternalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "IFstream.H"
#include "DynamicList.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::operationMode
>
Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::operationModeNames
({
    { operationMode::fixedPower, "power" },
    { operationMode::fixedHeatFlux, "flux" },
    { operationMode::fixedHeatTransferCoeff, "coefficient" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::readWaveRecord()
{
    // Try multiple paths to find the wave file
    fileName resolvedPath = waveFile_;
    
    if (!isFile(resolvedPath))
    {
        // Try relative to case directory
        resolvedPath = db().time().path()/waveFile_;
    }
    
    if (!isFile(resolvedPath))
    {
        FatalErrorInFunction
            << "waveMode = 'file' but cannot find waveFile for patch "
            << this->patch().name() << nl
            << "Tried paths:" << nl
            << "  " << waveFile_ << nl
            << "  " << db().time().path()/waveFile_ << nl
            << "Please set 'waveFile' to a valid path."
            << exit(FatalError);
    }

    Info<< "Reading wave record from: " << resolvedPath << endl;

    IFstream in(resolvedPath);

    if (!in.good())
    {
        FatalErrorInFunction
            << "Cannot open waveFile '" << resolvedPath << "' for patch "
            << this->patch().name() << nl
            << "Check path and permissions."
            << exit(FatalError);
    }

    DynamicList<scalar> tList, etaList;
    string line;

    while (in.getLine(line))
    {
        // Skip empty lines
        if (line.empty()) continue;
        
        // Trim leading whitespace
        size_t firstNonSpace = line.find_first_not_of(" \t");
        if (firstNonSpace == std::string::npos) continue;
        
        // Skip comment lines
        if (line[firstNonSpace] == '#' || line[firstNonSpace] == '%') continue;
        
        // Replace commas with spaces for CSV compatibility
        std::replace(line.begin(), line.end(), ',', ' ');
        
        // Parse time and eta
        IStringStream iss(line);
        scalar tVal, etaVal;
        
        if (iss >> tVal >> etaVal)
        {
            tList.append(tVal);
            etaList.append(etaVal);
        }
    }

    if (tList.size() < 2)
    {
        FatalErrorInFunction
            << "Wave record '" << resolvedPath << "' has fewer than 2 entries."
            << nl << "Need at least 2 (t, eta) pairs for interpolation."
            << exit(FatalError);
    }

    tSeries_.transfer(tList);
    etaSeries_.transfer(etaList);

    tMin_ = min(tSeries_);
    tMax_ = max(tSeries_);
    
    // Statistics
    scalar minEta = min(etaSeries_);
    scalar maxEta = max(etaSeries_);
    
    Info<< "  Loaded " << tSeries_.size() << " time points" << nl
        << "  Time range: [" << tMin_ << ", " << tMax_ << "] s" << nl
        << "  Eta range: [" << minEta*1000 << ", " << maxEta*1000 << "] mm" << endl;
}


Foam::scalar Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::etaAtTime
(
    const scalar t
) const
{
    if (waveMode_ == "none")
    {
        return 0.0;
    }

    if (waveMode_ == "monochromatic")
    {
        return 0.5*H_*Foam::cos(omega_*t + phase_);
    }

    if (waveMode_ == "file")
    {
        if (tSeries_.empty()) return 0.0;

        if (t <= tMin_) return etaSeries_[0];
        if (t >= tMax_) return etaSeries_.last();

        // Binary search for efficiency with large files
        label lo = 0;
        label hi = tSeries_.size() - 1;
        
        while (hi - lo > 1)
        {
            label mid = (lo + hi) / 2;
            if (tSeries_[mid] <= t)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }

        // Linear interpolation
        const scalar t0 = tSeries_[lo];
        const scalar t1 = tSeries_[hi];
        const scalar e0 = etaSeries_[lo];
        const scalar e1 = etaSeries_[hi];

        const scalar w = (t - t0)/(t1 - t0 + SMALL);
        return e0 + w*(e1 - e0);
    }

    WarningInFunction
        << "Unknown waveMode '" << waveMode_
        << "' -> assuming eta(t) = 0." << endl;

    return 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::
waveExternalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), temperatureCoupledBase::KMethodType::mtLookup),
    mode_(fixedHeatTransferCoeff),
    Q_(0),
    q_(p.size(), 0),
    h_(p.size(), 0),
    Ta_(),
    relaxation_(1),
    emissivity_(0),
    qrPrevious_(p.size(), 0),
    qrRelaxation_(1),
    qrName_("none"),
    thicknessLayers_(),
    kappaLayers_(),
    // Wave extensions
    stillWaterLevel_(0.20),
    hWater_(),   // autoPtr, initialized empty
    TaWater_(),  // autoPtr, initialized empty
    hAir_(10.0),
    TaAir_(260.0),
    waveMode_("none"),
    waveFile_(""),
    tSeries_(),
    etaSeries_(),
    tMin_(0.0),
    tMax_(0.0),
    H_(0.0),
    T_(1.0),
    omega_(2.0*constant::mathematical::pi),
    phase_(0.0)
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::
waveExternalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(operationModeNames.get("mode", dict)),
    Q_(0),
    q_(),
    h_(),
    Ta_(),
    relaxation_(dict.getOrDefault<scalar>("relaxation", 1)),
    emissivity_(dict.getOrDefault<scalar>("emissivity", 0)),
    qrPrevious_(),
    qrRelaxation_(dict.getOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    thicknessLayers_(),
    kappaLayers_(),
    // Wave extensions
    stillWaterLevel_(dict.getOrDefault<scalar>("stillWaterLevel", 0.20)),
    hWater_(Function1<scalar>::NewIfPresent("hWater", dict, word::null, &db())),
    TaWater_(Function1<scalar>::NewIfPresent("TaWater", dict, word::null, &db())),
    hAir_(dict.getOrDefault<scalar>("hAir", 10.0)),
    TaAir_(dict.getOrDefault<scalar>("TaAir", 260.0)),
    waveMode_(dict.getOrDefault<word>("waveMode", "none")),
    waveFile_(dict.getOrDefault<fileName>("waveFile", "constant/waveRecord.dat")),
    tSeries_(),
    etaSeries_(),
    tMin_(0.0),
    tMax_(0.0),
    H_(dict.getOrDefault<scalar>("H", 0.0)),
    T_(dict.getOrDefault<scalar>("T", 1.0)),
    omega_(2.0*constant::mathematical::pi / dict.getOrDefault<scalar>("T", 1.0)),
    phase_(dict.getOrDefault<scalar>("phase", 0.0))
{
    // Mode-specific initialization
    switch (mode_)
    {
        case fixedPower:
            Q_ = dict.get<scalar>("Q");
            break;

        case fixedHeatFlux:
            q_ = scalarField("q", dict, p.size());
            break;

        case fixedHeatTransferCoeff:
            h_ = scalarField("h", dict, p.size());
            Ta_ = Function1<scalar>::New("Ta", dict, &db());
            if (dict.found("thicknessLayers"))
            {
                dict.readEntry("thicknessLayers", thicknessLayers_);
                dict.readEntry("kappaLayers", kappaLayers_);
            }
            break;
    }

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (qrName_ != "none")
    {
        if (dict.found("qrPrevious"))
        {
            qrPrevious_ = scalarField("qrPrevious", dict, p.size());
        }
        else
        {
            qrPrevious_.setSize(p.size(), 0);
        }
    }

    if (dict.found("refValue"))
    {
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }

    // Read wave record if using file mode
    if (waveMode_ == "file")
    {
        readWaveRecord();
    }
    
    // Info output
    if (waveMode_ != "none")
    {
        Info<< "waveExternalWallHeatFluxTemperature on patch " << p.name() << ":" << nl
            << "  waveMode = " << waveMode_ << nl
            << "  stillWaterLevel = " << stillWaterLevel_ << " m" << nl
            << "  hWater = " << (hWater_ ? "Function1" : "not set")
            << ", TaWater = " << (TaWater_ ? "Function1" : "not set") << nl
            << "  hAir = " << hAir_ << " W/m2/K, TaAir = " << TaAir_ << " K" << endl;

        if (waveMode_ == "monochromatic")
        {
            Info<< "  H = " << H_ << " m, T = " << T_ << " s" << endl;
        }
    }
}


Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::
waveExternalWallHeatFluxTemperatureFvPatchScalarField
(
    const waveExternalWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mode_(ptf.mode_),
    Q_(ptf.Q_),
    q_(ptf.q_, mapper),
    h_(ptf.h_, mapper),
    Ta_(ptf.Ta_.clone()),
    relaxation_(ptf.relaxation_),
    emissivity_(ptf.emissivity_),
    qrPrevious_(ptf.qrPrevious_, mapper),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    // Wave extensions
    stillWaterLevel_(ptf.stillWaterLevel_),
    hWater_(ptf.hWater_.clone()),    // Clone the Function1
    TaWater_(ptf.TaWater_.clone()),  // Clone the Function1
    hAir_(ptf.hAir_),
    TaAir_(ptf.TaAir_),
    waveMode_(ptf.waveMode_),
    waveFile_(ptf.waveFile_),
    tSeries_(ptf.tSeries_),
    etaSeries_(ptf.etaSeries_),
    tMin_(ptf.tMin_),
    tMax_(ptf.tMax_),
    H_(ptf.H_),
    T_(ptf.T_),
    omega_(ptf.omega_),
    phase_(ptf.phase_)
{}


Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::
waveExternalWallHeatFluxTemperatureFvPatchScalarField
(
    const waveExternalWallHeatFluxTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    mode_(tppsf.mode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_.clone()),
    relaxation_(tppsf.relaxation_),
    emissivity_(tppsf.emissivity_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_),
    // Wave extensions
    stillWaterLevel_(tppsf.stillWaterLevel_),
    hWater_(tppsf.hWater_.clone()),   // Clone the Function1
    TaWater_(tppsf.TaWater_.clone()), // Clone the Function1
    hAir_(tppsf.hAir_),
    TaAir_(tppsf.TaAir_),
    waveMode_(tppsf.waveMode_),
    waveFile_(tppsf.waveFile_),
    tSeries_(tppsf.tSeries_),
    etaSeries_(tppsf.etaSeries_),
    tMin_(tppsf.tMin_),
    tMax_(tppsf.tMax_),
    H_(tppsf.H_),
    T_(tppsf.T_),
    omega_(tppsf.omega_),
    phase_(tppsf.phase_)
{}


Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::
waveExternalWallHeatFluxTemperatureFvPatchScalarField
(
    const waveExternalWallHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf),
    mode_(tppsf.mode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_.clone()),
    relaxation_(tppsf.relaxation_),
    emissivity_(tppsf.emissivity_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_),
    // Wave extensions
    stillWaterLevel_(tppsf.stillWaterLevel_),
    hWater_(tppsf.hWater_.clone()),   // Clone the Function1
    TaWater_(tppsf.TaWater_.clone()), // Clone the Function1
    hAir_(tppsf.hAir_),
    TaAir_(tppsf.TaAir_),
    waveMode_(tppsf.waveMode_),
    waveFile_(tppsf.waveFile_),
    tSeries_(tppsf.tSeries_),
    etaSeries_(tppsf.etaSeries_),
    tMin_(tppsf.tMin_),
    tMax_(tppsf.tMax_),
    H_(tppsf.H_),
    T_(tppsf.T_),
    omega_(tppsf.omega_),
    phase_(tppsf.phase_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.autoMap(mapper);
            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.autoMap(mapper);
            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.autoMap(mapper);
    }
}


void Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const waveExternalWallHeatFluxTemperatureFvPatchScalarField& ewhftpsf =
        refCast<const waveExternalWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.rmap(ewhftpsf.q_, addr);
            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.rmap(ewhftpsf.h_, addr);
            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.rmap(ewhftpsf.qrPrevious_, addr);
    }
}


void Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tp = *this;
    
    // Get current time for Function1 evaluation
    const scalar t = this->db().time().value();

    // Handle radiative flux
    scalarField qr(Tp.size(), Zero);

    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField>(qrName_);

        if (qrRelaxation_ < 1)
        {
            qr = qrRelaxation_*qr + (1 - qrRelaxation_)*qrPrevious_;
            qrPrevious_ = qr;
        }
    }

    // Store relaxation terms
    const scalarField valueFraction0 = valueFraction();
    const scalarField refValue0 = refValue();

    switch (mode_)
    {
        case fixedPower:
        {
            refGrad() = (Q_/gSum(patch().magSf()) + qr)/kappa(Tp);
            refValue() = Tp;
            valueFraction() = 0;
            break;
        }
        case fixedHeatFlux:
        {
            refGrad() = (q_ + qr)/kappa(Tp);
            refValue() = Tp;
            valueFraction() = 0;
            break;
        }
        case fixedHeatTransferCoeff:
        {
            // ================================================================
            // WAVE-DRIVEN PER-FACE h/Ta LOGIC
            // ================================================================

            // Current time and water level
            const scalar eta = etaAtTime(t);
            const scalar zWaterLevel = stillWaterLevel_ + eta;

            const vectorField& Cf = patch().Cf();

            // Effective per-face h and Ta
            scalarField hEff(h_.size());
            scalarField TaEff(h_.size());

            // Get hWater and TaWater values at current time (or defaults)
            const scalar hWaterValue = hWater_ ? hWater_->value(t) : 400.0;
            const scalar TaWaterValue = TaWater_ ? TaWater_->value(t) : 278.15;

            label nSubmerged = 0;

            forAll(Cf, i)
            {
                const scalar z = Cf[i].z();

                if (z <= zWaterLevel)
                {
                    // Submerged - use water properties
                    hEff[i] = (hWaterValue > 0) ? hWaterValue : h_[i];
                    TaEff[i] = TaWaterValue;
                    nSubmerged++;
                }
                else
                {
                    // Exposed to air
                    hEff[i] = (hAir_ > 0) ? hAir_ : h_[i];
                    TaEff[i] = TaAir_;
                }
            }

            // Debug output at output times
            if (this->db().time().outputTime())
            {
                const scalar fracSubmerged = (Cf.size() > 0)
                    ? scalar(nSubmerged)/scalar(Cf.size())*100.0
                    : 0.0;

                Info<< "waveExternalWallHeatFlux [" << patch().name()
                    << "] t=" << t << " s:" << nl
                    << "  eta = " << eta*1000 << " mm, "
                    << "water level = " << zWaterLevel*1000 << " mm" << nl
                    << "  hWater = " << hWaterValue << " W/m2/K, "
                    << "TaWater = " << TaWaterValue << " K" << nl
                    << "  submerged: " << nSubmerged << "/" << Cf.size()
                    << " (" << fracSubmerged << "%)" << endl;
            }

            // ================================================================
            // ORIGINAL THIN-LAYER LOGIC (using hEff/TaEff)
            // ================================================================

            scalar totalSolidRes = 0;
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    const scalar l = thicknessLayers_[iLayer];
                    if (kappaLayers_[iLayer] > 0)
                    {
                        totalSolidRes += l/kappaLayers_[iLayer];
                    }
                }
            }

            scalarField hp(hEff.size());
            forAll(hEff, i)
            {
                hp[i] = 1.0/(1.0/hEff[i] + totalSolidRes);
            }

            scalarField hpTa(hEff.size());
            forAll(hEff, i)
            {
                hpTa[i] = hp[i]*TaEff[i];
            }

            // Emissivity corrections (if enabled)
            if (emissivity_ > 0)
            {
                if (totalSolidRes > 0)
                {
                    scalarField TpLambda(hEff.size());
                    forAll(hEff, i)
                    {
                        TpLambda[i] = hEff[i]/(hEff[i] + 1.0/totalSolidRes);
                    }

                    scalarField Ts(hEff.size());
                    scalarField lambdaTa4(hEff.size());

                    forAll(hEff, i)
                    {
                        Ts[i] = TpLambda[i]*Tp[i]
                              + (1 - TpLambda[i])*TaEff[i];
                        lambdaTa4[i] = pow4((1 - TpLambda[i])*TaEff[i]);
                    }

                    forAll(hEff, i)
                    {
                        hp[i] += emissivity_*sigma.value()
                               *(pow4(Ts[i]) - lambdaTa4[i])/Tp[i];
                        hpTa[i] += emissivity_*sigma.value()
                                 *(lambdaTa4[i] + pow4(TaEff[i]));
                    }
                }
                else
                {
                    forAll(hEff, i)
                    {
                        hp[i] += emissivity_*sigma.value()*pow3(Tp[i]);
                        hpTa[i] += emissivity_*sigma.value()*pow4(TaEff[i]);
                    }
                }
            }

            const scalarField kappaDeltaCoeffs
            (
                this->kappa(Tp)*patch().deltaCoeffs()
            );

            refGrad() = 0;

            forAll(Tp, i)
            {
                if (qr[i] < 0)
                {
                    const scalar hpmqr = hp[i] - qr[i]/Tp[i];

                    refValue()[i] = hpTa[i]/hpmqr;
                    valueFraction()[i] = hpmqr/(hpmqr + kappaDeltaCoeffs[i]);
                }
                else
                {
                    refValue()[i] = (hpTa[i] + qr[i])/hp[i];
                    valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
                }
            }

            break;
        }
    }

    // Apply relaxation
    valueFraction() =
        relaxation_*valueFraction()
      + (1 - relaxation_)*valueFraction0;

    refValue() = relaxation_*refValue() + (1 - relaxation_)*refValue0;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("mode", operationModeNames[mode_]);
    temperatureCoupledBase::write(os);

    switch (mode_)
    {
        case fixedPower:
        {
            os.writeEntry("Q", Q_);
            break;
        }
        case fixedHeatFlux:
        {
            q_.writeEntry("q", os);
            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.writeEntry("h", os);
            if (Ta_)
            {
                Ta_->writeData(os);
            }

            if (relaxation_ < 1)
            {
                os.writeEntry("relaxation", relaxation_);
            }

            if (emissivity_ > 0)
            {
                os.writeEntry("emissivity", emissivity_);
            }

            if (thicknessLayers_.size())
            {
                thicknessLayers_.writeEntry("thicknessLayers", os);
                kappaLayers_.writeEntry("kappaLayers", os);
            }

            break;
        }
    }

    os.writeEntry("qr", qrName_);

    if (qrName_ != "none")
    {
        os.writeEntry("qrRelaxation", qrRelaxation_);
        qrPrevious_.writeEntry("qrPrevious", os);
    }

    // Wave extension parameters
    os.writeEntry("waveMode", waveMode_);
    os.writeEntry("stillWaterLevel", stillWaterLevel_);

    // Write hWater - supports both scalar and Function1
    if (hWater_)
    {
        hWater_->writeData(os);
    }

    // Write TaWater - supports both scalar and Function1
    if (TaWater_)
    {
        TaWater_->writeData(os);
    }

    os.writeEntry("hAir", hAir_);
    os.writeEntry("TaAir", TaAir_);
    
    if (waveMode_ == "file")
    {
        os.writeEntry("waveFile", waveFile_);
    }
    else if (waveMode_ == "monochromatic")
    {
        os.writeEntry("H", H_);
        os.writeEntry("T", T_);
        os.writeEntry("phase", phase_);
    }

    refValue().writeEntry("refValue", os);
    refGrad().writeEntry("refGradient", os);
    valueFraction().writeEntry("valueFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        waveExternalWallHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //