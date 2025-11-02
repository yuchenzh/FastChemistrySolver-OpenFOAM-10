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

\*---------------------------------------------------------------------------*/

#include "FastChemistryModel.H"
#include "janafThermo.H"
#include "perfectGas.H"
#include "specie.H"
#include "sutherlandTransport.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Define the single thermo type used by FastChemistryModel
// Using the same structure as typedefThermo macro
typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > FastChemistryModelThermoType;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline ChemistryModelHelper::ChemistryModelHelper
(
    const basicChemistryModel& chemistry
)
:
    fastChemPtr_(const_cast<void*>(static_cast<const void*>(&chemistry)))
{
    IOdictionary chemDict
    (
        IOobject
        (
            chemistry.thermo().phasePropertyName("chemistryProperties"),
            chemistry.thermo().phasePropertyName(),
            chemistry.thermo().T().mesh().time().constant(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    Switch fastODE = chemDict.lookupOrDefault<Switch> ("optimizedODE", false);
    if (fastODE)
    {
        if (! chemDict.found("fastODESolver"))
        {
            FatalErrorInFunction
                << "optimizedODE is true but fastODESolver not found in chemistryProperties"
                << exit(FatalError);
        }
        else
        {
            word solverName = chemDict.lookup<word>("fastODESolver");
        
            const word chemSolverName = solverName + "<" + "FastChemistryModel" + 
                "<" + chemistry.thermo().thermoName() + ">>";
            Info << "Selecting fast chemistry solver " << chemSolverName << endl;
        }



    }

}


inline scalarField ChemistryModelHelper::getRRGivenYTP
(
    const scalarField& Y,
    const scalar T,
    const scalar p,
    const scalar deltaT,
    scalar& deltaTChem,
    const scalar& rho,
    const scalar& rho0
) const
{
    // Direct cast to FastChemistryModel (only one thermo type supported)
    const FastChemistryModel<FastChemistryModelThermoType>* fastChem =
        static_cast<const FastChemistryModel<FastChemistryModelThermoType>*>(fastChemPtr_);

    return fastChem->getRRGivenYTP(Y, T, p, deltaT, deltaTChem, rho, rho0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
