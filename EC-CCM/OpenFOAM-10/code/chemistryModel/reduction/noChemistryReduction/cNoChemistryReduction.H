/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::CCMchemistryReductionMethods::none

Description
    A CCMchemistryReductionMethod which does nothing to allow reduction to be
    switched-off.

SourceFiles
    cNoChemistryReduction.C

\*---------------------------------------------------------------------------*/

#ifndef cNoChemistryReduction_H
#define cNoChemistryReduction_H

#include "CCMchemistryReductionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CCMchemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                            Class none Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class none
:
    public CCMchemistryReductionMethod<ThermoType>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        none
        (
            const IOdictionary& dict,
            CCMchemistryModel<ThermoType>& chemistry
        );


    //- Destructor
    virtual ~none();


    // Member Functions

        //- Return mechanism reduction inactive
        virtual bool active() const
        {
            return false;
        }

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            List<label>& ctos,
            DynamicList<label>& stoc,
            const label li
        );

        virtual void update()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CCMchemistryReductionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cNoChemistryReduction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
