/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::CCMchemistryReductionMethod

Description
    An abstract class for methods of chemical mechanism reduction

SourceFiles
    CCMchemistryReductionMethod.C
    CCMchemistryReductionMethods.C
    CCMchemistryReductionMethodNew.C
    CCMchemistryReductionMethodI.H

\*---------------------------------------------------------------------------*/

#ifndef CCMchemistryReductionMethod_H
#define CCMchemistryReductionMethod_H

#include "IOdictionary.H"
#include "DynamicField.H"
#include "Switch.H"
#include "cpuTime.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ThermoType>
class CCMchemistryModel;

/*---------------------------------------------------------------------------*\
                   Class CCMchemistryReductionMethod Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class CCMchemistryReductionMethod
{
protected:

    // Protected Data

        //- Dictionary that store the algorithm data
        const dictionary coeffsDict_;

        //- Reference to the chemistry model
        CCMchemistryModel<ThermoType>& chemistry_;

        //- Total number of species
        const label nSpecie_;

        //- Number of active species
        label nActiveSpecies_;

        //- List of disabled reactions (disabled = true)
        Field<bool> reactionsDisabled_;

        //- List of active species (active = true)
        List<bool> activeSpecies_;


    //- Protected Member Functions

        //- Initialise reduction of the mechanism
        void initReduceMechanism();

        //- End reduction of the mechanism
        void endReduceMechanism(List<label>& ctos, DynamicList<label>& stoc);


private:

    // Private Data

        //- Switch to select performance logging
        Switch log_;

        //- Tolerance for the mechanism reduction algorithm
        scalar tolerance_;

        //- CPU time for logging
        cpuTime cpuTime_;

        //- ...
        int64_t sumnActiveSpecies_;

        //- ...
        int64_t sumn_;

        //- ...
        scalar reduceMechCpuTime_;

        // Log file for the average time spent reducing the chemistry
        autoPtr<OFstream> cpuReduceFile_;

        // Write average number of species
        autoPtr<OFstream> nActiveSpeciesFile_;


public:

    //- Runtime type information
    TypeName("CCMchemistryReductionMethod");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CCMchemistryReductionMethod,
        dictionary,
        (
            const IOdictionary& dict,
            CCMchemistryModel<ThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        CCMchemistryReductionMethod
        (
            CCMchemistryModel<ThermoType>& chemistry
        );

        //- Construct from components
        CCMchemistryReductionMethod
        (
            const IOdictionary& dict,
            CCMchemistryModel<ThermoType>& chemistry
        );


    // Selector

        static autoPtr<CCMchemistryReductionMethod<ThermoType>> New
        (
            const IOdictionary& dict,
            CCMchemistryModel<ThermoType>& chemistry
        );


    //- Destructor
    virtual ~CCMchemistryReductionMethod();


    // Member Functions

        //- Return mechanism reduction active
        virtual bool active() const
        {
            return true;
        }

        //- Return the number of species
        inline label nSpecie();

        //- Return the number of active species
        inline label nActiveSpecies() const;

        //- Return the tolerance
        inline scalar tolerance() const;

        //- Return the active species
        inline const List<bool>& activeSpecies() const;

        //- Return whether or not a reaction is disabled
        inline bool reactionDisabled(const label i) const;

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            List<label>& ctos,
            DynamicList<label>& stoc,
            const label li
        ) = 0;

        //- ...
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CCMchemistryReductionMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CCMchemistryReductionMethod.C"
    #include "CCMchemistryReductionMethodNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
