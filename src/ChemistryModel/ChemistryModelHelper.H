/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::ChemistryModelHelper

Description
    Helper class to access FastChemistryModel methods without modifying
    the base class interface.

    FastChemistryModel only supports one thermo type:
    sutherlandTransport + sensibleEnthalpy + janafThermo + perfectGas + specie

    Usage example:
    \verbatim
        autoPtr<basicChemistryModel> chemistry =
            basicChemistryModel::New(thermo);

        ChemistryModelHelper helper(chemistry());

        scalarField RR = helper.getRRGivenYTP
        (
            Y, T, p, dt, deltaTChem, rho, rho
        );
    \endverbatim

SourceFiles
    ChemistryModelHelperI.H (inline functions)

\*---------------------------------------------------------------------------*/

#ifndef ChemistryModelHelper_H
#define ChemistryModelHelper_H

#include "basicChemistryModel.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
template<class ThermoType> class FastChemistryModel;

/*---------------------------------------------------------------------------*\
                    Class ChemistryModelHelper Declaration
\*---------------------------------------------------------------------------*/

class ChemistryModelHelper
{
    // Private data

        //- Pointer to FastChemistryModel (cached after first cast)
        mutable void* fastChemPtr_;


public:

    // Constructors

        //- Construct from chemistry model reference
        inline ChemistryModelHelper(const basicChemistryModel& chemistry);


    // Member Functions

        //- Get reaction rates given Y, T, P
        //  Returns the reaction rate for each species [kg/m^3/s]
        inline scalarField getRRGivenYTP
        (
            const scalarField& Y,
            const scalar T,
            const scalar p,
            const scalar deltaT,
            scalar& deltaTChem,
            const scalar& rho,
            const scalar& rho0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ChemistryModelHelperI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
