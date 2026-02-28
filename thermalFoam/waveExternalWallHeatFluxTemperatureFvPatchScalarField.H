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

Class
    Foam::waveExternalWallHeatFluxTemperatureFvPatchScalarField

Description
    Mixed thermal boundary condition with wave-driven per-face switching
    between water and air convection properties.

    Wave modes:
      - none: No wave forcing (original behavior)
      - file: Read eta(t) from external file
      - monochromatic: Synthetic regular wave eta = H/2 * cos(omega*t + phase)

Usage
    \verbatim
    bluffFace
    {
        type            waveExternalWallHeatFluxTemperature;
        mode            coefficient;
        
        Ta              constant 278.15;
        h               uniform 100;
        
        waveMode        file;
        waveFile        "constant/waveEta.dat";
        stillWaterLevel 0.18;
        
        hWater          350;      // or use tableFile for time-varying
        TaWater         278.15;   // or use tableFile for time-varying
        hAir            12;
        TaAir           268.15;

        kappaMethod     lookup;
        kappa           Kth;

        value           uniform 263.15;
    }
    \endverbatim

    hWater and TaWater can be specified as:
      - A constant scalar: hWater 350;
      - A Function1 (e.g., tableFile):
        hWater
        {
            type            tableFile;
            file            "constant/water_htc.txt";
            outOfBounds     clamp;
            interpolationScheme linear;
        }

SourceFiles
    waveExternalWallHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef waveExternalWallHeatFluxTemperatureFvPatchScalarField_H
#define waveExternalWallHeatFluxTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "Function1.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class waveExternalWallHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class waveExternalWallHeatFluxTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    // Public data

        //- Operation mode enumeration
        enum operationMode
        {
            fixedPower,
            fixedHeatFlux,
            fixedHeatTransferCoeff
        };

        //- Operation mode names
        static const Enum<operationMode> operationModeNames;


private:

    // Private Data

        //- Operation mode
        operationMode mode_;

        //- Heat power [W]
        scalar Q_;

        //- Heat flux [W/m²]
        scalarField q_;

        //- Heat transfer coefficient [W/m²K]
        scalarField h_;

        //- Ambient temperature [K]
        autoPtr<Function1<scalar>> Ta_;

        //- Relaxation for the wall temperature (thermal inertia)
        scalar relaxation_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivity_;

        //- Cache qr for relaxation
        scalarField qrPrevious_;

        //- Relaxation for qr
        scalar qrRelaxation_;

        //- Name of the radiative heat flux
        const word qrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        // ============================================================
        // Wave-aware extensions (Omonigbehin, INRS-ETE)
        // ============================================================

        //- Still water level [m]
        scalar stillWaterLevel_;

        //- Water-side heat transfer coefficient [W/m²K] - supports Function1
        autoPtr<Function1<scalar>> hWater_;

        //- Water-side ambient temperature [K] - supports Function1
        autoPtr<Function1<scalar>> TaWater_;

        //- Air-side heat transfer coefficient [W/m²K]
        scalar hAir_;

        //- Air-side ambient temperature [K]
        scalar TaAir_;

        //- Wave mode: "none", "file", "monochromatic"
        word waveMode_;

        //- Wave record file path
        fileName waveFile_;

        //- Time series from file [s]
        scalarField tSeries_;

        //- Eta series from file [m]
        scalarField etaSeries_;

        //- Minimum time in record
        scalar tMin_;

        //- Maximum time in record
        scalar tMax_;

        //- Wave height for monochromatic mode [m]
        scalar H_;

        //- Wave period for monochromatic mode [s]
        scalar T_;

        //- Angular frequency [rad/s]
        scalar omega_;

        //- Initial phase [rad]
        scalar phase_;


    // Private Member Functions

        //- Read wave record from file
        void readWaveRecord();

        //- Interpolate eta at given time
        scalar etaAtTime(const scalar t) const;


public:

    //- Runtime type information
    TypeName("waveExternalWallHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        waveExternalWallHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveExternalWallHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        waveExternalWallHeatFluxTemperatureFvPatchScalarField
        (
            const waveExternalWallHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        waveExternalWallHeatFluxTemperatureFvPatchScalarField
        (
            const waveExternalWallHeatFluxTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new waveExternalWallHeatFluxTemperatureFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        waveExternalWallHeatFluxTemperatureFvPatchScalarField
        (
            const waveExternalWallHeatFluxTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new waveExternalWallHeatFluxTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //