/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::CombustionHelpers

Description
    Chemistry-specific calculation functions for CCM chemistry model including
    element mass fractions, equivalence ratio, and scalar dissipation rate.

SourceFiles
    CombustionHelpers.C

\*---------------------------------------------------------------------------*/

#ifndef CombustionHelpers_H
#define CombustionHelpers_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "Switch.H"
#include "scalarList.H"
#include "hashedWordList.H"
#include "fluidReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CombustionHelpers Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class CombustionHelpers
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to thermo
        const fluidReactionThermo& thermo_;

        //- Reference to CCM dictionary
        const dictionary& CCMdict_;

        //- Reference to mass fraction fields
        const PtrList<volScalarField>& Yvf_;

        //- Reference to update variables list
        const hashedWordList& updateVars_;

        //- References to J coefficient arrays
        const scalarList& JhCoeff_;
        const scalarList& JcCoeff_;
        const scalarList& JnCoeff_;
        const scalarList& JoCoeff_;
        const scalarList& Jh_h2oCoeff_;
        const scalarList& Jc_co2Coeff_;
        const scalarList& Jo_co2_h2oCoeff_;

        //- References to element mass fraction fields
        volScalarField& Jh_;
        volScalarField& Jc_;
        volScalarField& Jn_;
        volScalarField& Jo_;
        volScalarField& Jh_h2o_;
        volScalarField& Jc_co2_;
        volScalarField& Jo_co2_h2o_;
        volScalarField& J_;

        //- References to equivalence ratio and scalar dissipation rate fields
        volScalarField& phieq_;
        volScalarField& chi_;

        //- Reference to style switches
        const Switch& oldStyleChi_;
        const Switch& oldStylePhi_;

        //- Reference to fuel oxygen to carbon ratio
        const scalar& fuelOtoC_;

public:

    // Constructors

        //- Construct from references to required data
        CombustionHelpers
        (
            const fvMesh& mesh,
            const fluidReactionThermo& thermo,
            const dictionary& CCMdict,
            const PtrList<volScalarField>& Yvf,
            const hashedWordList& updateVars,
            const scalarList& JhCoeff,
            const scalarList& JcCoeff,
            const scalarList& JnCoeff,
            const scalarList& JoCoeff,
            const scalarList& Jh_h2oCoeff,
            const scalarList& Jc_co2Coeff,
            const scalarList& Jo_co2_h2oCoeff,
            volScalarField& Jh,
            volScalarField& Jc,
            volScalarField& Jn,
            volScalarField& Jo,
            volScalarField& Jh_h2o,
            volScalarField& Jc_co2,
            volScalarField& Jo_co2_h2o,
            volScalarField& J,
            volScalarField& phieq,
            volScalarField& chi,
            const Switch& oldStyleChi,
            const Switch& oldStylePhi,
            const scalar& fuelOtoC
        );

        //- Disallow default bitwise copy construction
        CombustionHelpers(const CombustionHelpers&) = delete;


    //- Destructor
    ~CombustionHelpers() = default;


    // Member Functions

        //- Update carbon element mass fraction
        void updateJc();

        //- Update hydrogen element mass fraction
        void updateJh();

        //- Update oxygen element mass fraction
        void updateJo();

        //- Update nitrogen element mass fraction
        void updateJn();

        //- Update hydrogen element mass fraction excluding H2O
        void updateJh_h2o();

        //- Update carbon element mass fraction excluding CO2
        void updateJc_co2();

        //- Update oxygen element mass fraction excluding CO2 and H2O
        void updateJo_co2_h2o();

        //- Update equivalence ratio
        void updatePhi();

        //- Update scalar dissipation rate
        void updateChi();

        //- Update principal variables (excluding species mass fractions and temperature)
        void updatePV();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CombustionHelpers&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CombustionHelpers.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //