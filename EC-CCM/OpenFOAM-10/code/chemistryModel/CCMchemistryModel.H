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
    Foam::CCMchemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source terms
    with optional support for TDAC mechanism reduction and tabulation.

    References:
    \verbatim
        Contino, F., Jeanmart, H., Lucchini, T., & D’Errico, G. (2011).
        Coupling of in situ adaptive tabulation and dynamic adaptive chemistry:
        An effective method for solving combustion in engine simulations.
        Proceedings of the Combustion Institute, 33(2), 3057-3064.

        Contino, F., Lucchini, T., D'Errico, G., Duynslaegher, C.,
        Dias, V., & Jeanmart, H. (2012).
        Simulations of advanced combustion modes using detailed chemistry
        combined with tabulation and mechanism reduction techniques.
        SAE International Journal of Engines,
        5(2012-01-0145), 185-196.

        Contino, F., Foucher, F., Dagaut, P., Lucchini, T., D’Errico, G., &
        Mounaïm-Rousselle, C. (2013).
        Experimental and numerical analysis of nitric oxide effect on the
        ignition of iso-octane in a single cylinder HCCI engine.
        Combustion and Flame, 160(8), 1476-1483.

        Contino, F., Masurier, J. B., Foucher, F., Lucchini, T., D’Errico, G., &
        Dagaut, P. (2014).
        CFD simulations using the TDAC method to model iso-octane combustion
        for a large range of ozone seeding and temperature conditions
        in a single cylinder HCCI engine.
        Fuel, 137, 179-184.
    \endverbatim

SourceFiles
    CCMchemistryModelI.H
    CCMchemistryModel.C

\*---------------------------------------------------------------------------*/
#define  CCMINFO Info << "#====== CCM Information ======#" << endl;
#define  STOP throw "error";
#define SHUTDOWN FatalErrorIn(__FUNCTION__) << " Shut Down By the User To Debug" << abort(FatalError);

// Timing macros for CCM performance analysis with parallel synchronization
#define CCM_SYNC() \
    if (Pstream::parRun()) { scalar sync = 1.0; reduce(sync, sumOp<scalar>()); }

#define CCM_TIMING_START(stepName, debugTimeFlag, timer) \
    if (debugTimeFlag) \
    { \
        CCM_SYNC(); \
        timer.cpuTimeIncrement(); \
    }

#define CCM_TIMING_END(stepName, debugTimeFlag, timer, chemistryModel) \
    if (debugTimeFlag) \
    { \
        CCM_SYNC(); \
        label stepIdx = (chemistryModel).getNextStepIndex(); \
        stepNames_.append(#stepName); \
        (chemistryModel).stepTimes_[stepIdx] = timer.cpuTimeIncrement(); \
    }

#ifndef CCMchemistryModel_H
#define CCMchemistryModel_H

#include "CCModeChemistryModel.H"
#include "ReactionList.H"
#include "ODESystem.H"
#include "volFields.H"
#include "multiComponentMixture.H"
#include "CCMchemistryReductionMethod.H"
#include "CCMchemistryTabulationMethod.H"
#include "DynamicField.H"
#include "CCMUtilities.H"
#include "CCMDebug.H"
#include "ParallelComm.H"
#include "CombustionHelpers.H"
#include "atomicWeights.H"
#include "ReactionEntry.H"
#include "RREntry.H"
#include "RoundRobin.H"
#include "basicFastChemistryModel.H"

// Helper classes
#include "CCMHelpers.H"

// Hash table type definitions for pre-allocation optimization
using rrEntries = HashTable<RREntry,string>;
using rrEntriesList = List<rrEntries>;
#include "IOmanip.H"
#include "cpuTime.H"
#include "Pstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CCMchemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class CCMchemistryModel
:
    public CCModeChemistryModel
{
    // Private classes

        //- Class to define scope of reaction evaluation. Runs pre-evaluate
        //  hook on all reactions on construction and post-evaluate on
        //  destruction.
        class reactionEvaluationScope
        {
            const CCMchemistryModel<ThermoType>& chemistry_;

        public:

            reactionEvaluationScope
            (
                const CCMchemistryModel<ThermoType>& chemistry
            )
            :
                chemistry_(chemistry)
            {
                forAll(chemistry_.reactions_, i)
                {
                    chemistry_.reactions_[i].preEvaluate();
                }
            }

            ~reactionEvaluationScope()
            {
                forAll(chemistry_.reactions_, i)
                {
                    chemistry_.reactions_[i].postEvaluate();
                }
            }
        };


    // Private data

        //- Switch to select performance logging
        Switch log_;

        //- Switch to enable loadBalancing performance logging
        Switch loadBalancing_;

        //- Type of the Jacobian to be calculated
        const jacobianType jacobianType_;

        //- Reference to the multi component mixture
        const multiComponentMixture<ThermoType>& mixture_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermos_;

        //- Reactions
        const ReactionList<ThermoType> reactions_;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField::Internal> RR_;

        //- Temporary mass fraction field
        mutable scalarField Y_;

        //- Temporary simplified mechanism mass fraction field
        DynamicField<scalar> sY_;

        //- Temporary concentration field
        mutable scalarField c_;

        //- Temporary simplified mechanism concentration field
        DynamicField<scalar> sc_;

        //- Specie-temperature-pressure workspace fields
        mutable FixedList<scalarField, 5> YTpWork_;

        //- Specie-temperature-pressure workspace matrices
        mutable FixedList<scalarSquareMatrix, 2> YTpYTpWork_;

        //- Mechanism reduction method
        autoPtr<CCMchemistryReductionMethod<ThermoType>> mechRedPtr_;

        //- Mechanism reduction method reference
        CCMchemistryReductionMethod<ThermoType>& mechRed_;

        //- Tabulation method
        autoPtr<CCMchemistryTabulationMethod> tabulationPtr_;

        //- Tabulation method reference
        CCMchemistryTabulationMethod& tabulation_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;

        //- The ratio of oxygen to carbon (mole) for the fuel
        scalar fuelOtoC_;

        //- Jh coefficient for each specie in the reaction mechanism
        scalarList JhCoeff_;

        //- Jc coefficient for each specie in the reaction mechanism
        scalarList JcCoeff_;

        //- Jn coefficient for each specie in the reaction mechanism
        scalarList JnCoeff_;

        //- Jo coefficient for each specie in the reaction mechanism
        scalarList JoCoeff_;

        //- Jh_h2o Jh coefficient for each specie excluding H2O
        scalarList Jh_h2oCoeff_;

        //- Jc_co2 Jc coefficient for each specie excluding CO2
        scalarList Jc_co2Coeff_;

        //- Jo_co2_h2o Jo coefficient for each specie excluding CO2 and H2O
        scalarList Jo_co2_h2oCoeff_;

        //- CCM dictionary
        dictionary CCMdict_;

        //- CCM Utilities helper class
        mutable autoPtr<CCMUtilities<ThermoType>> ccmUtilities_;

        //- CCM Debug helper class
        mutable autoPtr<CCMDebug<ThermoType>> ccmDebug_;

        //- CCM Parallel Communication helper class
        mutable autoPtr<ParallelComm<ThermoType>> parallelComm_;

        //- CCM Combustion helper class
        mutable autoPtr<CombustionHelpers<ThermoType>> combustionHelpers_;

        //- A list of pricipal variables
        hashedWordList principalVars_;

        //- A list of species-only principal variables (subset of principalVars_)
        hashedWordList speciesPrincipalVars_;

        //- A list of regular variables
        hashedWordList   regularVars_;

        //- A list of principal variables to be updated
        hashedWordList updateVars_;

        //- Error control variables (subset of principalVars for encoding)
        hashedWordList ecVars_;


        //- The number of slices between min and max values for most 
        //-  regular variables
        label nSlice_;

        //- The threshold to ignore for a min value
        scalar ignoreMin_;

        //- Min/max/span of the principal variables
        HashTable<scalar, word> pvMin_;
        HashTable<scalar, word> pvMax_;
        HashTable<scalar, word> pvSpan_;

        //- Several element mass fraction fields (Jh, Jc, Jn, Jo ...)
        volScalarField Jh_,Jc_,Jn_,Jo_,Jh_h2o_,Jc_co2_,Jo_co2_h2o_, J_;

        //- Equivalence ratio(phi) and scalar dissipation rate(chi)
        volScalarField phieq_, chi_;

        //- The maximum indices representations
        label maxiRepresentation_;

        //- RoundRobin Communication
        RoundRobin RoundRobinCommunicator_;

        //- Old style scalar dissipation rate
        Switch oldStyleChi_;

        //- Old style equivalence ratio
        Switch oldStylePhi_;

        // //- Acceleration Information output
        // OFstream accelerationInfo_;
        
        //- Flag to set debug mode
        Switch debugMode_;

        //- Flag to examine Ydiff details
        Switch examineYdiff_;

        //- Flag to shutdown immediately
        Switch shutdownImmediately_;

        //- Flag to use optimized communication pattern
        Switch optimizedCommunication_;

        //- Ydiff for debug
        PtrList<volScalarField> Ydiff_;

        //- Temporary global maximum mass fraction field
        scalarField Ymax_;

        //- Flag to enable detailed timing analysis
        Switch debugTime_;

        //- reaction index (named zoneIndex in the program)
        List<string> zoneIndex_;

        //- remaidner of the reaction index (named zoneRemainder in the program)
        labelList zoneRemainder_;

        //- reaction entries for each core to calculate
        reactionEntries gatheredReactionEntries_;

        //- Timing arrays for the solve function steps
        mutable scalarField stepTimes_;

        //- Step names corresponding to timing entries
        mutable wordList stepNames_;

        //- CPU timer for step timing
        mutable cpuTime stepTimer_;

        //- Current step index for automatic timing
        mutable label currentStepIndex_;

        //- Grouping Hashtable time
        mutable cpuTime groupingTimer_;

        //- Grouping time scalar
        scalar groupingTime_;

        //- EC mode settings
        Switch ecEnabled_;
        //wordList ecInitialVars_;
        label numECVarsToAdd_;
        label numECVarsToRemove_;
        label ecUpdateFreq_;
        label currentStep_;
        bool ecInitialized_;

        //- Flag to treat pressure variations in high Mach cases
        Switch highMach_;

        //- Pre-allocated hash tables for performance optimization
        //  These are sized in constructor and cleared before each use
        
        //- Pre-allocated hash tables for cell grouping (Steps 2-3)
        mutable reactionEntriesList preAllocatedToCore_;
        mutable reactionEntriesList preAllocatedFromCore_;
        
        //- Pre-allocated hash table for local reaction rate storage (Step 4) 
        mutable rrEntries preAllocatedLocalRREntries_;
        
        //- Pre-allocated hash tables for optimized communication (Step 6)
        mutable rrEntriesList preAllocatedReturnToCore_;
        mutable rrEntriesList preAllocatedReceivedRREntries_;

        //- Optimized ODE solver
        const autoPtr<basicFastChemistryModel> fastChemistryPtr_;

    // Private Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField::Internal>& RR();



        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("CCM");

    //- Structure to return variable ranking information
    struct VariableRanking 
    {
        scalarField importance;     // Normalized importance for each principalVar species
        labelList candidateIndices; // PrincipalVar indices that can be added to ecVars
        labelList sortedIndices;    // Species indices sorted by importance (descending)
        labelList sortedECIndices;  // Species indices sorted (only for ecVars_)
        bool underControl;         // True if optimization is necessary (candidates exist)
        
        VariableRanking() : underControl(false) {}
    };
    // Constructors

        //- Construct from thermo
        CCMchemistryModel(const fluidReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        CCMchemistryModel(const CCMchemistryModel&) = delete;


    //- Destructor
    virtual ~CCMchemistryModel();


    // Member Functions

        //- Return reference to the mixture
        inline const multiComponentMixture<ThermoType>& mixture() const;

        //- The reactions
        inline const PtrList<Reaction<ThermoType>>& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermos() const;

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m^3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );

            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;


        // ODE functions (overriding abstract functions in ODE.H)

            virtual void derivatives
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt,
                scalarSquareMatrix& J
            ) const;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& Y,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;

            // Functions to distribute reactionEntry/solve/update the reaction rates
            VariableRanking distributeReactionEntry(
                const hashedWordList& encodingVars = hashedWordList(),
                bool returnRanking = false
            );

            //- Find n least important variables from current ecVars_
            //  Returns list of variable names to remove
            hashedWordList findLeastImportantEcVars(const VariableRanking& ranking, label n) const;

            void updateReactionRate();

            void distributeReactionRate();

            void initPvMinMaxSpan();

            //- Calculate the reaction rates give Y,T and P. 
            // Note, density should be provided!
            // deltaTChem will be modified!
            scalarField getRRGivenYTP
            (
                scalarField& Y,
                scalar& T, 
                scalar& p,
                const scalar& deltaT,
                scalar& deltaTChem,
                const scalar& rho,
                const scalar& rho0
            );


            // //- Reset timing step counter
            // void resetTimingSteps() const;

            //- Get next step index and increment counter
            label getNextStepIndex() const;

            //- update zoneIndex_ and zoneRemainder_ if mesh has changed
            void updateCCM4MeshChange();


            // //- Add step name to the dynamic list
            // void addStepName(const word& stepName) const;


        //- Calculate variable ranking based on std analysis
        VariableRanking calculateVariableRanking(const hashedWordList& currentVars, scalar tolerance);

        // Mechanism reduction access functions

            //- Return true if specie i is active
            inline bool active(const label i) const;

            //- Set specie i active
            inline void setActive(const label i);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CCMchemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CCMchemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CCMchemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
