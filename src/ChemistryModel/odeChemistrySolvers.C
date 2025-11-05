/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "OptRosenbrock34.H"
#include "OptRodas34.H"
#include "OptSeulex.H"
#include "FastChemistryModel.H"
#include "basicFastChemistryModel.H"
// #include "forGases.H"
// #include "forLiquids.H"
//#include "makeChemistrySolver.H"
// #include "forThermo.H"
//#include "chemistryModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// #define defineMyChemistrySolver(Transport,Energy,Thermo,Equation,Specie,solver)  \
//                                                                                \
//     typedefThermo(Transport,Energy,Thermo,Equation,Specie);                    \
//                                                                                \
//                                                                                \
//     defineChemistrySolver                                                      \
//     (                                                                          \
//         FastChemistryModel,                                                      \
//         Transport##Energy##Thermo##Equation##Specie                            \
//     );                                                                         \
//     makeChemistrySolver                                                        \
//     (                                                                          \
//         solver,                                                                 \
//         FastChemistryModel,                                                      \
//         Transport##Energy##Thermo##Equation##Specie                            \
//     );                                                                         \
//                                                                                \


//template instantiation

// Define FastChemistryModel type name (non-template class)
defineTypeNameAndDebug(FastChemistryModel, 0);

// Simplified macro for non-template FastChemistryModel
#define makeChemistrySolver(Solver)                                            \
    typedef Solver<FastChemistryModel> Solver##FastChemistryModel;             \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Solver##FastChemistryModel,                                            \
        (                                                                      \
            word(Solver##FastChemistryModel::typeName_())                      \
          + "<FastChemistryModel>"                                             \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        basicFastChemistryModel,                                               \
        Solver##FastChemistryModel,                                            \
        mesh                                                                 \
    )

// Register ODE solvers
makeChemistrySolver(OptRosenbrock34);
makeChemistrySolver(OptRodas34);
makeChemistrySolver(OptSeulex);






//typedefThermo(sutherlandTransport,   sensibleEnthalpy,   janafThermo,    perfectGas, specie);
//defineChemistrySolver(chemistryModel, sutherlandTransportsensibleEnthalpyjanafThermoperfectGasspecie);
//makeChemistrySolver(OptRosenbrock34,chemistryModel,sutherlandTransportsensibleEnthalpyjanafThermoperfectGasspecie);
}


// ************************************************************************* //
