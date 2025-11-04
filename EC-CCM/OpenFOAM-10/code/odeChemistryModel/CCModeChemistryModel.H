/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::CCModeChemistryModel

Description
    Extends base chemistry model adding an ODESystem and the reduction maps
    needed for tabulation.

SourceFiles
    CCModeChemistryModelI.H
    CCModeChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef CCModeChemistryModel_H
#define CCModeChemistryModel_H

#include "basicChemistryModel.H"
#include "ODESystem.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ThermoType>
class CCMchemistryModel;

template<class ThermoType>
class CCMchemistryReductionMethod;


/*---------------------------------------------------------------------------*\
                     Class CCModeChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class CCModeChemistryModel
:
    public basicChemistryModel,
    public ODESystem
{
    // Private data

        //- Reference to the field of specie mass fractions
        const PtrList<volScalarField>& Yvf_;

        //- Number of species
        label nSpecie_;

        //- Is chemistry reduction active
        bool reduction_;

        //- Temporary map from complete to simplified concentration fields
        //  c -> sc
        List<label> cTos_;

        //- Temporary map from simplified to complete concentration fields
        //  sc -> c
        DynamicList<label> sToc_;


public:

    template<class ThermoType>
    friend class CCMchemistryModel;

    template<class ThermoType>
    friend class CCMchemistryReductionMethod;

    //- Runtime type information
    TypeName("CCModeChemistryModel");


    // Constructors

        //- Construct from thermo
        CCModeChemistryModel(const fluidReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        CCModeChemistryModel(const CCModeChemistryModel&) = delete;


    //- Destructor
    virtual ~CCModeChemistryModel();


    // Member Functions

        //- Create and return a TDAC log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

        //- The number of species
        inline virtual label nSpecie() const;

        //- Allow the reduction method to reset the number of species
        inline void setNSpecie(const label newNs);

        //- Number of ODE's to solve
        inline virtual label nEqns() const;

        //- Return a reference to the list of mass fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return true if mechanism reduction is active
        inline bool reduction() const;

        //- Return the index in the complete set of species
        //  corresponding to the index si in the simplified set of species
        inline label sToc(const label si) const;

        //- Return the index in the simplified set of species
        //  corresponding to the index si in the complete set of species
        inline label cTos(const label ci) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CCModeChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CCModeChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
