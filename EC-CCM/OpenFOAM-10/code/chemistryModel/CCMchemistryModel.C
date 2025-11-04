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

\*---------------------------------------------------------------------------*/

#include "CCMchemistryModel.H"
#include "UniformField.H"
#include "localEulerDdtScheme.H"
#include "cpuLoad.H"
#include "cNoChemistryTabulation.H"
#include <iomanip>
#include "cNoChemistryReduction.H"
#include <vector>


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::CCMchemistryModel<ThermoType>::CCMchemistryModel
(
    const fluidReactionThermo& thermo
)
:
    CCModeChemistryModel(thermo),
    log_(this->lookupOrDefault("log", false)),
    loadBalancing_(this->lookupOrDefault("loadBalancing", false)),
    jacobianType_
    (
        this->found("jacobian")
      ? jacobianTypeNames_.read(this->lookup("jacobian"))
      : jacobianType::fast
    ),
    mixture_(refCast<const multiComponentMixture<ThermoType>>(this->thermo())),
    specieThermos_(mixture_.specieThermos()),
    reactions_(mixture_.species(), specieThermos_, this->mesh(), *this),
    RR_(nSpecie_),
    Y_(nSpecie_),
    c_(nSpecie_),
    YTpWork_(scalarField(nSpecie_ + 2)),
    YTpYTpWork_(scalarSquareMatrix(nSpecie_ + 2)),
    mechRedPtr_
    (
        autoPtr<CCMchemistryReductionMethod<ThermoType>>
        (
            new CCMchemistryReductionMethods::none<ThermoType>
            (
                *this,
                *this
            )
        )
    ),
    mechRed_(*mechRedPtr_),
    tabulationPtr_
    (
        autoPtr<CCMchemistryTabulationMethod>
        (
            new CCMchemistryTabulationMethods::none
            (
                *this,
                *this
            )
            
        )
    ),
    tabulation_(*tabulationPtr_),
    fuelOtoC_(readScalar(this->subDict("CCM").lookup("ratioOxygenToCarbonElementInFuel"))),
    JhCoeff_(scalarList(nSpecie_, 0.)),
    JcCoeff_(scalarList(nSpecie_, 0.)),
    JnCoeff_(scalarList(nSpecie_, 0.)),
    JoCoeff_(scalarList(nSpecie_, 0.)),
    Jh_h2oCoeff_(scalarList(nSpecie_, 0.)),
    Jc_co2Coeff_(scalarList(nSpecie_, 0.)),
    Jo_co2_h2oCoeff_(scalarList(nSpecie_, 0.)),
    CCMdict_(this->subDict("CCM")),
    ccmUtilities_(nullptr),
    ccmDebug_(nullptr),
    parallelComm_(nullptr),
    combustionHelpers_(nullptr),
    principalVars_(CCMdict_.lookup<hashedWordList>("principalVars")),
    speciesPrincipalVars_(wordList()),
    regularVars_(wordList()),
    updateVars_(wordList()),
    nSlice_(CCMdict_.lookupOrDefault<label>("nSlice",50)),
    ignoreMin_(CCMdict_.lookupOrDefault<scalar>("ignoreMin",1e-6)),
    Jh_
    (
        volScalarField
        (
            IOobject
            (
                "Jh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jc_
    (
        volScalarField
        (
            IOobject
            (
                "Jc",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jn_
    (
        volScalarField
        (
            IOobject
            (
                "Jn",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jo_
    (
        volScalarField
        (
            IOobject
            (
                "Jo",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jh_h2o_
    (
        volScalarField
        (
            IOobject
            (
                "Jh_h2o",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jc_co2_
    (
        volScalarField
        (
            IOobject
            (
                "Jc_co2",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jo_co2_h2o_
    (
        volScalarField
        (
            IOobject
            (
                "Jo_co2_h2o",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    J_
    (
        volScalarField
        (
            IOobject
            (
                "J",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    phieq_
    (
        volScalarField
        (
            IOobject
            (
                "phieq",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    chi_
    (
        volScalarField
        (
            IOobject
            (
                "chi",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimTime, 0)
        )
    ),
    maxiRepresentation_(29728),
    RoundRobinCommunicator_(),
    oldStyleChi_(CCMdict_.lookupOrDefault<Switch>("oldStyleChi",false)),
    oldStylePhi_(CCMdict_.lookupOrDefault<Switch>("oldStylePhi",false)),
    // accelerationInfo_("CCMAccInfo"),
    debugMode_(CCMdict_.lookupOrDefault<Switch>("debugMode", false)),
    examineYdiff_(CCMdict_.lookupOrDefault<Switch>("examineYdiff", false)),
    shutdownImmediately_(CCMdict_.lookupOrDefault<Switch>("shutdownImmediately", false)),
    optimizedCommunication_(CCMdict_.lookupOrDefault<Switch>("optimizedCommunication", true)),
    Ydiff_(nSpecie_),
    Ymax_(nSpecie_),
    debugTime_(CCMdict_.lookupOrDefault<Switch>("debugTime", false)),
    zoneIndex_(mesh().nCells()),
    zoneRemainder_(mesh().nCells()),
    gatheredReactionEntries_(0),
    stepTimes_(16, 0.0),
    stepNames_(),
    stepTimer_(),
    currentStepIndex_(0),
    groupingTimer_(),
    groupingTime_(0.0),
    ecEnabled_(CCMdict_.found("ecMode") ? Switch(CCMdict_.subDict("ecMode").lookupOrDefault<bool>("enabled", false)) : Switch(false)),
    numECVarsToAdd_(CCMdict_.found("ecMode") ? CCMdict_.subDict("ecMode").lookupOrDefault<label>("numECVarsToAdd", 3) : 3),
    numECVarsToRemove_(CCMdict_.found("ecMode") ? CCMdict_.subDict("ecMode").lookupOrDefault<label>("numECVarsToRemove", 3) : 3),
    ecUpdateFreq_(CCMdict_.found("ecMode") ? CCMdict_.subDict("ecMode").lookupOrDefault<label>("updateFreq", 10) : 10),
    currentStep_(0),
    ecInitialized_(false),
    highMach_(CCMdict_.lookupOrDefault<Switch>("highMach", false)),
    preAllocatedToCore_(),
    preAllocatedFromCore_(),
    preAllocatedLocalRREntries_(),
    preAllocatedReturnToCore_(),
    preAllocatedReceivedRREntries_(),
    fastChemistryPtr_(
        basicFastChemistryModel::New(thermo)
    )
{

    
    // Pre-allocate hash tables with double the principalVars size for efficiency
    pvMin_.resize(2 * principalVars_.size());
    pvMax_.resize(2 * principalVars_.size());
    pvSpan_.resize(2 * principalVars_.size());
    
    // Initialize speciesPrincipalVars_ and ecInitialVars_.
    // ecVars_ is initialized with any principalVars_ that are not species
    forAll(principalVars_, pi)
    {
        const word& varName = principalVars_[pi];
        if (mixture_.species().found(varName))
        {
            speciesPrincipalVars_.append(varName);
        }
        else
        {
            ecVars_.append(varName);
        }
    }
    ecInitialized_ = true;

    if (debugMode_)
    {
        Info << "Enabling debug mode" << endl;
    }    

    if (examineYdiff_)
    {
        Info << "Output Ydiff examination" << endl;
    }

    if (shutdownImmediately_)
    {
        Info << "Shutdown immediately after outputing Ydiff" << endl;
    }

    if (oldStyleChi_)
    {
        Info << "Using old style scalar dissipation rate" << endl;
    }
    else
    {
        Info << "Using correct scalar dissipation rate" << endl;
    }

    if (oldStylePhi_)
    {
        Info << "Using the equivalence ratio that will NOT change during the combustion process" << endl;
    }
    else
    {
        Info << "Using the equivalence ratio that will change during the combustion process" << endl;
    }

    if (optimizedCommunication_)
    {
        Info << "Using optimized targeted communication pattern" << endl;
    }
    else
    {
        Info << "Using original broadcast communication pattern" << endl;
    }

    if (highMach_)
    {
        Info << "Treating pressure variations in high Mach cases" << endl;
    }
    


    // accelerationInfo_.precision(5);
    // accelerationInfo_ << setw(8) << "time" << tab 
    //     << setw(8) << "totalCells" << tab 
    //     << setw(8) << "phaseCells" << tab 
    //     << setw(8) << "accRatio" << endl;

    // Initialize CCM Utilities helper
    ccmUtilities_.set
    (
        new CCMUtilities<ThermoType>
        (
            mixture_,
            specieThermos_,
            CCMdict_,
            nSpecie_,
            principalVars_,
            updateVars_,
            regularVars_,
            oldStylePhi_,
            JhCoeff_,
            JcCoeff_,
            JnCoeff_,
            JoCoeff_,
            Jh_h2oCoeff_,
            Jc_co2Coeff_,
            Jo_co2_h2oCoeff_,
            Yvf_,
            Ymax_,
            pvMin_,
            pvMax_,
            pvSpan_,
            ignoreMin_,
            this->mesh()
        )
    );

    // Initialize CCM Debug helper
    ccmDebug_.set
    (
        new CCMDebug<ThermoType>
        (
            nSpecie_,
            Ymax_,
            mixture_,
            stepTimes_,
            stepNames_,
            currentStepIndex_
        )
    );

    // Initialize CCM Parallel Communication helper
    parallelComm_.set
    (
        new ParallelComm<ThermoType>
        (
            RoundRobinCommunicator_
        )
    );

    // Initialize CCM Combustion helper
    combustionHelpers_.set
    (
        new CombustionHelpers<ThermoType>
        (
            this->mesh(),
            this->thermo(),
            CCMdict_,
            Yvf_,
            updateVars_,
            JhCoeff_,
            JcCoeff_,
            JnCoeff_,
            JoCoeff_,
            Jh_h2oCoeff_,
            Jc_co2Coeff_,
            Jo_co2_h2oCoeff_,
            Jh_,
            Jc_,
            Jn_,
            Jo_,
            Jh_h2o_,
            Jc_co2_,
            Jo_co2_h2o_,
            J_,
            phieq_,
            chi_,
            oldStyleChi_,
            oldStylePhi_,
            fuelOtoC_
        )
    );

    // Initialize Element mass fractions (Jc, Jh, Jn, Jo ...)
    ccmUtilities_->initJCoeffs();


    // Initialize the regular variables, that min/max/span are required
    ccmUtilities_->initRegularVars();

    // Initialize variable lists, including the list of principal variables (PV) 
    // and the PVs to be updated
    ccmUtilities_->initPvLists();

    // initialize min/max/span for PVs
    ccmUtilities_->initPvMinMaxSpan();


    // initalized Ydiff if debug
    if (debugMode_ && examineYdiff_)
    {
        forAll(Ydiff_, i)
        {
            Ydiff_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Ydiff." + Yvf_[i].name(),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar(dimless, 0)
                )
            );
        }
    }

    CCMINFO;

    // examine min/max for PVs


    Info << nl;
    Info << "Principal variables: ";
    ccmUtilities_->outputHashedWordList(principalVars_);
    // forAll(principalVars_, pi){Info << principalVars_[pi] << tab;} Info << endl;

    Info << "Regular variables: ";
    ccmUtilities_->outputHashedWordList(regularVars_);
    // forAll(regularVars_, ri){Info << regularVars_[ri] << tab;} Info << endl;

    Info << "Variables to be update: ";
    ccmUtilities_->outputHashedWordList(updateVars_);
    // forAll(updateVars_, ui){Info << updateVars_[ui] << tab;} Info << endl;

    Info << "nSlice: " << nSlice_ << endl;
    Info << "ignoreMin: " << ignoreMin_ << endl;

    Info << "Error control variables(ecVars_): ";
    ccmUtilities_->outputHashedWordList(ecVars_);

    combustionHelpers_->updatePV();

    ccmUtilities_->examineMinMax();

    combustionHelpers_->updatePV();

    // For all principal variables, set the write option to AUTO_WRITE
    forAll(principalVars_, pi)
    {
        volScalarField& field = this->mesh().objectRegistry::lookupObjectRef<volScalarField>(principalVars_[pi]);
        field.writeOpt() = IOobject::AUTO_WRITE;
    }





    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Yvf_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.T().mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "CCMchemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction() << endl;

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialised (by default 'active' is true)
    if (reduction_)
    {
        const basicSpecieMixture& composition = this->thermo().composition();

        forAll(Yvf_, i)
        {
            typeIOobject<volScalarField> header
            (
                Yvf_[i].name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ
            );

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
            if (!header.headerOk())
            {
                composition.setInactive(i);
            }
        }
    }

    // Create communicator
    if (Pstream::parRun())
    {
        string mode = CCMdict_.subDict("communicator").lookupOrDefault<word>("mode", "global");
        if (mode == "global")
        {
            Info << "Global communicator" << endl;
            RoundRobinCommunicator_ = RoundRobin(Pstream::nProcs());
        }
        else if (mode == "distributed")
        {
            Info << "Distributed communicator" << endl;
            label localCores = CCMdict_.subDict("communicator").lookupOrDefault<label>("localCores", 2);
            if (localCores < 2)
            {
                FatalErrorIn("CCMchemistryModel::CCMchemistryModel")
                << "Number of local cores must be greater than 1"
                << abort(FatalError);
            }
            
            Info << "localCores: " << localCores << endl;
            RoundRobinCommunicator_ = RoundRobin(Pstream::nProcs(), localCores);
        }
        
    }
    else // create a fake one for single core cases
    {
        RoundRobinCommunicator_ = RoundRobin(2);
    }
    RoundRobinCommunicator_.recordOpponents(false);

    // Pre-allocate hash tables for performance optimization
    // Size based on mesh and processor configuration to prevent expansion
    label contactNum = (Pstream::parRun()) ? RoundRobinCommunicator_.nCores_ : 1;
    label estimatedCellsPerCore = mesh().nCells() / contactNum;
    
    // Pre-allocate toCore and fromCore lists with sufficient capacity
    // Use 2x estimated size to handle load imbalance and prevent expansion
    label hashTableCapacity = max(1024, 2 * estimatedCellsPerCore);
    
    preAllocatedToCore_.resize(contactNum);
    preAllocatedFromCore_.resize(contactNum);
    preAllocatedReturnToCore_.resize(contactNum);
    preAllocatedReceivedRREntries_.resize(contactNum);
    
    forAll(preAllocatedToCore_, i)
    {
        preAllocatedToCore_[i] = reactionEntries(hashTableCapacity);
        preAllocatedFromCore_[i] = reactionEntries(hashTableCapacity);
        preAllocatedReturnToCore_[i] = rrEntries(hashTableCapacity);
        preAllocatedReceivedRREntries_[i] = rrEntries(hashTableCapacity);
    }
    
    // Pre-allocate local RR entries with estimated phase space size
    // Use larger capacity to handle worst-case phase space expansion
    label localRRCapacity = max(2048, 4 * estimatedCellsPerCore);
    preAllocatedLocalRREntries_ = rrEntries(localRRCapacity);
    
    Info << "Hash tables pre-allocated with capacity: " << hashTableCapacity 
         << " (contact cores: " << contactNum << ", local RR: " << localRRCapacity << ")" << endl;

    if (log_)
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }
  
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::CCMchemistryModel<ThermoType>::~CCMchemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::CCMchemistryModel<ThermoType>::derivatives
(
    const scalar time,
    const scalarField& YTp,
    const label li,
    scalarField& dYTpdt
) const
{
    if (reduction_)
    {
        forAll(sToc_, i)
        {
            Y_[sToc_[i]] = max(YTp[i], 0);
        }
    }
    else
    {
        forAll(Y_, i)
        {
            Y_[i] = max(YTp[i], 0);
        }
    }

    const scalar T = YTp[nSpecie_];
    const scalar p = YTp[nSpecie_ + 1];

    // Evaluate the mixture density
    scalar rhoM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        rhoM += Y_[i]/specieThermos_[i].rho(p, T);
    }
    rhoM = 1/rhoM;

    // Evaluate the concentrations
    for (label i=0; i<Y_.size(); i ++)
    {
        c_[i] = rhoM/specieThermos_[i].W()*Y_[i];
    }

    // Evaluate contributions from reactions
    dYTpdt = Zero;
    forAll(reactions_, ri)
    {
        if (!mechRed_.reactionDisabled(ri))
        {
            reactions_[ri].dNdtByV
            (
                p,
                T,
                c_,
                li,
                dYTpdt,
                reduction_,
                cTos_,
                0
            );
        }
    }

    // Reactions return dNdtByV, so we need to convert the result to dYdt
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar WiByrhoM = specieThermos_[sToc(i)].W()/rhoM;
        scalar& dYidt = dYTpdt[i];
        dYidt *= WiByrhoM;
    }

    // Evaluate the effect on the thermodynamic system ...

    // Evaluate the mixture Cp
    scalar CpM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        CpM += Y_[i]*specieThermos_[i].Cp(p, T);
    }

    // dT/dt
    scalar& dTdt = dYTpdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        dTdt -= dYTpdt[i]*specieThermos_[sToc(i)].Ha(p, T);
    }
    dTdt /= CpM;

    // dp/dt = 0 (pressure is assumed constant)
    scalar& dpdt = dYTpdt[nSpecie_ + 1];
    dpdt = 0;
}


template<class ThermoType>
void Foam::CCMchemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& YTp,
    const label li,
    scalarField& dYTpdt,
    scalarSquareMatrix& J
) const
{
    if (reduction_)
    {
        forAll(sToc_, i)
        {
            Y_[sToc_[i]] = max(YTp[i], 0);
        }
    }
    else
    {
        forAll(c_, i)
        {
            Y_[i] = max(YTp[i], 0);
        }
    }

    const scalar T = YTp[nSpecie_];
    const scalar p = YTp[nSpecie_ + 1];

    // Evaluate the specific volumes and mixture density
    scalarField& v = YTpWork_[0];
    for (label i=0; i<Y_.size(); i++)
    {
        v[i] = 1/specieThermos_[i].rho(p, T);
    }
    scalar rhoM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        rhoM += Y_[i]*v[i];
    }
    rhoM = 1/rhoM;

    // Evaluate the concentrations
    for (label i=0; i<Y_.size(); i ++)
    {
        c_[i] = rhoM/specieThermos_[i].W()*Y_[i];
    }

    // Evaluate the derivatives of concentration w.r.t. mass fraction
    scalarSquareMatrix& dcdY = YTpYTpWork_[0];
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar rhoMByWi = rhoM/specieThermos_[sToc(i)].W();
        switch (jacobianType_)
        {
            case jacobianType::fast:
                {
                    dcdY(i, i) = rhoMByWi;
                }
                break;
            case jacobianType::exact:
                for (label j=0; j<nSpecie_; j++)
                {
                    dcdY(i, j) =
                        rhoMByWi*((i == j) - rhoM*v[sToc(j)]*Y_[sToc(i)]);
                }
                break;
        }
    }

    // Evaluate the mixture thermal expansion coefficient
    scalar alphavM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        alphavM += Y_[i]*rhoM*v[i]*specieThermos_[i].alphav(p, T);
    }

    // Evaluate contributions from reactions
    dYTpdt = Zero;
    scalarSquareMatrix& ddNdtByVdcTp = YTpYTpWork_[1];
    for (label i=0; i<nSpecie_ + 2; i++)
    {
        for (label j=0; j<nSpecie_ + 2; j++)
        {
            ddNdtByVdcTp[i][j] = 0;
        }
    }
    forAll(reactions_, ri)
    {
        if (!mechRed_.reactionDisabled(ri))
        {
            reactions_[ri].ddNdtByVdcTp
            (
                p,
                T,
                c_,
                li,
                dYTpdt,
                ddNdtByVdcTp,
                reduction_,
                cTos_,
                0,
                nSpecie_,
                YTpWork_[1],
                YTpWork_[2]
            );
        }
    }

    // Reactions return dNdtByV, so we need to convert the result to dYdt
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar WiByrhoM = specieThermos_[sToc(i)].W()/rhoM;
        scalar& dYidt = dYTpdt[i];
        dYidt *= WiByrhoM;

        for (label j=0; j<nSpecie_; j++)
        {
            scalar ddNidtByVdYj = 0;
            switch (jacobianType_)
            {
                case jacobianType::fast:
                    {
                        const scalar ddNidtByVdcj = ddNdtByVdcTp(i, j);
                        ddNidtByVdYj = ddNidtByVdcj*dcdY(j, j);
                    }
                    break;
                case jacobianType::exact:
                    for (label k=0; k<nSpecie_; k++)
                    {
                        const scalar ddNidtByVdck = ddNdtByVdcTp(i, k);
                        ddNidtByVdYj += ddNidtByVdck*dcdY(k, j);
                    }
                    break;
            }

            scalar& ddYidtdYj = J(i, j);
            ddYidtdYj = WiByrhoM*ddNidtByVdYj + rhoM*v[sToc(j)]*dYidt;
        }

        scalar ddNidtByVdT = ddNdtByVdcTp(i, nSpecie_);
        for (label j=0; j<nSpecie_; j++)
        {
            const scalar ddNidtByVdcj = ddNdtByVdcTp(i, j);
            ddNidtByVdT -= ddNidtByVdcj*c_[sToc(j)]*alphavM;
        }

        scalar& ddYidtdT = J(i, nSpecie_);
        ddYidtdT = WiByrhoM*ddNidtByVdT + alphavM*dYidt;

        scalar& ddYidtdp = J(i, nSpecie_ + 1);
        ddYidtdp = 0;
    }

    // Evaluate the effect on the thermodynamic system ...

    // Evaluate the mixture Cp and its derivative
    scalarField& Cp = YTpWork_[3];
    scalar CpM = 0, dCpMdT = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        Cp[i] = specieThermos_[i].Cp(p, T);
        CpM += Y_[i]*Cp[i];
        dCpMdT += Y_[i]*specieThermos_[i].dCpdT(p, T);
    }

    // dT/dt
    scalarField& Ha = YTpWork_[4];
    scalar& dTdt = dYTpdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        Ha[sToc(i)] = specieThermos_[sToc(i)].Ha(p, T);
        dTdt -= dYTpdt[i]*Ha[sToc(i)];
    }
    dTdt /= CpM;

    // dp/dt = 0 (pressure is assumed constant)
    scalar& dpdt = dYTpdt[nSpecie_ + 1];
    dpdt = 0;

    // d(dTdt)/dY
    for (label i=0; i<nSpecie_; i++)
    {
        scalar& ddTdtdYi = J(nSpecie_, i);
        ddTdtdYi = 0;
        for (label j=0; j<nSpecie_; j++)
        {
            const scalar ddYjdtdYi = J(j, i);
            ddTdtdYi -= ddYjdtdYi*Ha[sToc(j)];
        }
        ddTdtdYi -= Cp[sToc(i)]*dTdt;
        ddTdtdYi /= CpM;
    }

    // d(dTdt)/dT
    scalar& ddTdtdT = J(nSpecie_, nSpecie_);
    ddTdtdT = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar dYidt = dYTpdt[i];
        const scalar ddYidtdT = J(i, nSpecie_);
        ddTdtdT -= dYidt*Cp[sToc(i)] + ddYidtdT*Ha[sToc(i)];
    }
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT /= CpM;

    // d(dTdt)/dp = 0 (pressure is assumed constant)
    scalar& ddTdtdp = J(nSpecie_, nSpecie_ + 1);
    ddTdtdp = 0;

    // d(dpdt)/dYiTp = 0 (pressure is assumed constant)
    for (label i=0; i<nSpecie_ + 2; i++)
    {
        scalar& ddpdtdYiTp = J(nSpecie_ + 1, i);
        ddpdtdYiTp = 0;
    }
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::CCMchemistryModel<ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    scalarField& tc = ttc.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    if (this->chemistry_)
    {
        reactionEvaluationScope scope(*this);

        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            const scalar Ti = T[celli];
            const scalar pi = p[celli];

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Yvf_[i][celli]/specieThermos_[i].W();
            }

            // A reaction's rate scale is calculated as it's molar
            // production rate divided by the total number of moles in the
            // system.
            //
            // The system rate scale is the average of the reactions' rate
            // scales weighted by the reactions' molar production rates. This
            // weighting ensures that dominant reactions provide the largest
            // contribution to the system rate scale.
            //
            // The system time scale is then the reciprocal of the system rate
            // scale.
            //
            // Contributions from forward and reverse reaction rates are
            // handled independently and identically so that reversible
            // reactions produce the same result as the equivalent pair of
            // irreversible reactions.

            scalar sumW = 0, sumWRateByCTot = 0;
            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];
                scalar omegaf, omegar;
                R.omega(pi, Ti, c_, celli, omegaf, omegar);

                scalar wf = 0;
                forAll(R.rhs(), s)
                {
                    wf += R.rhs()[s].stoichCoeff*omegaf;
                }
                sumW += wf;
                sumWRateByCTot += sqr(wf);

                scalar wr = 0;
                forAll(R.lhs(), s)
                {
                    wr += R.lhs()[s].stoichCoeff*omegar;
                }
                sumW += wr;
                sumWRateByCTot += sqr(wr);
            }

            tc[celli] =
                sumWRateByCTot == 0 ? vGreat : sumW/sumWRateByCTot*sum(c_);
        }
    }

    ttc.ref().correctBoundaryConditions();

    return ttc;
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::CCMchemistryModel<ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        reactionEvaluationScope scope(*this);

        scalarField& Qdot = tQdot.ref();

        forAll(Yvf_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = specieThermos_[i].Hf();
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::CCMchemistryModel<ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    tmp<volScalarField::Internal> tRR
    (
        volScalarField::Internal::New
        (
            "RR",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    reactionEvaluationScope scope(*this);

    scalar omegaf, omegar;

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Yvf_[i][celli];
            c_[i] = rhoi*Yi/specieThermos_[i].W();
        }

        const Reaction<ThermoType>& R = reactions_[ri];
        const scalar omegaI = R.omega(pi, Ti, c_, celli, omegaf, omegar);

        forAll(R.lhs(), s)
        {
            if (si == R.lhs()[s].index)
            {
                RR[celli] -= R.lhs()[s].stoichCoeff*omegaI;
            }
        }

        forAll(R.rhs(), s)
        {
            if (si == R.rhs()[s].index)
            {
                RR[celli] += R.rhs()[s].stoichCoeff*omegaI;
            }
        }

        RR[celli] *= specieThermos_[si].W();
    }

    return tRR;
}


template<class ThermoType>
void Foam::CCMchemistryModel<ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField& dNdtByV = YTpWork_[0];

    reactionEvaluationScope scope(*this);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Yvf_[i][celli];
            c_[i] = rhoi*Yi/specieThermos_[i].W();
        }

        dNdtByV = Zero;

        forAll(reactions_, ri)
        {
            if (!mechRed_.reactionDisabled(ri))
            {
                reactions_[ri].dNdtByV
                (
                    pi,
                    Ti,
                    c_,
                    celli,
                    dNdtByV,
                    reduction_,
                    cTos_,
                    0
                );
            }
        }

        for (label i=0; i<mechRed_.nActiveSpecies(); i++)
        {
            RR_[sToc(i)][celli] = dNdtByV[i]*specieThermos_[sToc(i)].W();
        }
    }
}


template<class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::CCMchemistryModel<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{   
    // forAll(thermo().T(),i)
    // {   
    //     const volScalarField& T = this->thermo().T();
    //     const volScalarField& p = this->thermo().p();
    //     const volScalarField& rho = this->thermo().rho();
    //     scalar dt = mesh().time().deltaTValue();
    //     scalar deltaTChem = 1e-10;
    //     label celli = 1;
    //     scalarField Y(nSpecie_);
    //     forAll(Y, i)
    //     {
    //         Y[i] = Yvf_[i][celli];
    //     }
    //     scalar Tvalue = T[celli];
    //     scalar pvalue = p[celli];
    //     scalar rhovalue = rho[celli];

    //     scalarField Yorg(Y); scalar Torg(Tvalue); scalar porg(pvalue); scalar rhoorg(rhovalue);
    //     scalarField RR = getRRGivenYTP(Y, Tvalue, pvalue, dt, deltaTChem, rhovalue,rhovalue);
        
    //     // check difference
    //     scalarField Ydiff = Yorg - Y; scalar Tdiff = Torg - Tvalue; scalar pdiff = porg - pvalue; scalar rhodiff = rhoorg - rhovalue;
    //     scalar threshold = 1e-10;
    //     if (Foam::max(mag(Ydiff)) > threshold || mag(Tdiff) > threshold || mag(pdiff) > threshold || mag(rhodiff) > threshold)
    //     {
    //         Info << "At cell " << celli << ", after first chemistry solve:" << endl;
    //         Info << "Ydiff max = " << Foam::max(mag(Ydiff)) << ", Tdiff = " << mag(Tdiff) << ", pdiff = " << mag(pdiff) << ", rhodiff = " << mag(rhodiff) << endl;
    //     }

    //     deltaTChem = 1e-10;
    //     scalarField RRFast = fastChemistryPtr_->getRRGivenYTP(Y, Tvalue, pvalue, dt, deltaTChem, rhovalue, rhovalue);
    //     scalarField RRDiff = RR - RRFast;
    //     forAll(RRDiff, i)
    //     {
    //         Info << "Species " << i 
    //                 << ": RR = " << RR[i] 
    //                 << ", RRFast = " << RRFast[i] 
    //                 << ", Diff = " << RRDiff[i] << endl;
    //     }
    //     throw "error";
    // }
    
    // throw "error";

    scalar deltaTMin = great;
    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    // update zoneIndex_ and zoneRemainder_ if mesh has changed
    updateCCM4MeshChange();



    // Reset step timing array
    if (debugTime_)
    {
        stepTimes_ = 0.0;

        // Reset automatic step counter
        currentStepIndex_ = 0;
        stepNames_.clear();

        // Reset timer
        stepTimer_.cpuTimeIncrement(); // Reset timer
    }


    
    if (ecEnabled_ && (currentStep_ % ecUpdateFreq_ == 0))
    {
        
        Info << "EC update at step " << currentStep_ << " with ecVars: ";
        ccmUtilities_->outputHashedWordList(ecVars_);

        VariableRanking ranking = this->distributeReactionEntry(ecVars_, true);
        bool ecOkay = ranking.underControl;
        
        if (!ecOkay)
        {
            while (!ecOkay)
            {
                label toAdd = min(numECVarsToAdd_, ranking.candidateIndices.size());
                Info << "The error is not under control" << endl;
                Info << "Adding " << toAdd << " variables:" << endl;
                
                for (label i = 0; i < toAdd; i++)
                {
                    word varName = speciesPrincipalVars_[ranking.candidateIndices[i]];
                    ecVars_.append(varName);
                    Info << "  " << varName << endl;
                }
                
                Info << "Updated ecVars: ";
                ccmUtilities_->outputHashedWordList(ecVars_);

                ranking = this->distributeReactionEntry(ecVars_, true);
                ecOkay = ranking.underControl;
            }
        }
        else
        {
            // try remove n least important variables
            Info << "Error has already been controlled" << endl;
            Info << "Trying to remove " << numECVarsToRemove_ << " least important variables:" << endl;
            hashedWordList varsToRemove = findLeastImportantEcVars(ranking, numECVarsToRemove_);
            
            // Report variables to be removed
            if (varsToRemove.size() > 0)
            {
                ccmUtilities_->outputHashedWordList(varsToRemove);
                
                hashedWordList newEcVars;
                forAll(ecVars_, i)
                {
                    if (! varsToRemove.found(ecVars_[i]))
                    {
                        newEcVars.append(ecVars_[i]);
                    }
                }
                ecVars_ = newEcVars;
                ranking = this->distributeReactionEntry(ecVars_, true);
                ecOkay = ranking.underControl;
            }
            else
            {
                Info << "  No variables to remove (safety limit reached)" << endl;
            }
            
            if (!ecOkay)
            {
                Info << "EC still not under control after removing variables, "
                     << "keeping current ecVars: " << endl;
                forAll(varsToRemove, vari)
                {
                    ecVars_.append(varsToRemove[vari]);
                }
                ranking = this->distributeReactionEntry(ecVars_, true);
                ecOkay = ranking.underControl;

                ccmUtilities_->outputHashedWordList(ecVars_);
            }
            else
            {
                Info << "EC under control after removing variables, "
                     << "final ecVars: ";
                ccmUtilities_->outputHashedWordList(ecVars_);
            }
        }
        
    }
    else if (ecEnabled_)
    {

        this->distributeReactionEntry(ecVars_, false);
    }
    else
    {
        // Normal mode - use principalVars
        this->distributeReactionEntry();
    }

    updateReactionRate();

    distributeReactionRate();

    deltaTMin = min(deltaTChem_).value();

    // Step 9: Cleanup & Statistics
    CCM_TIMING_START(Cleanup_Statistics, debugTime_, stepTimer_);
    

    mechRed_.update();
    tabulation_.update();
    
    const basicSpecieMixture& composition = this->thermo().composition();
    if (reduction_ && Pstream::parRun())
    {
        List<bool> active(composition.active());
        Pstream::listCombineGather(active, orEqOp<bool>());
        Pstream::listCombineScatter(active);

        forAll(active, i)
        {
            if (active[i])
            {
                composition.setActive(i);
            }
        }
    }
    CCM_TIMING_END(Cleanup_Statistics, debugTime_, stepTimer_, *this);

    // Output timing analysis if enabled
    if (debugTime_)
    {
        ccmDebug_->outputTimingAnalysis(stepTimes_);
    }

    // Increment step counter
    currentStep_++;

    return deltaTMin;
}


template<class ThermoType>
Foam::scalar Foam::CCMchemistryModel<ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ThermoType>
Foam::scalar Foam::CCMchemistryModel<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ThermoType>
typename Foam::CCMchemistryModel<ThermoType>::VariableRanking 
Foam::CCMchemistryModel<ThermoType>::calculateVariableRanking(
    const hashedWordList& currentVars, 
    scalar tolerance
)
{
    VariableRanking result;
    
    // Calculate importance scores (reuse debugStdRank logic)
    scalarField Ymaxstd(nSpecie_, 0.0);
    for (auto it = gatheredReactionEntries_.begin(); it != gatheredReactionEntries_.end(); it++)
    {
        const ReactionEntry& re = it();
        forAll(re.Y, yi)
        {
            if (re.count > 1)
            {
                Ymaxstd[yi] = max(Ymaxstd[yi], re.Ystd[yi]);
            }
        }
    }
    reduce(Ymaxstd, maxOp<scalarList>());
    
    // Normalize by Ymax
    forAll(Ymax_, yi)
    {
        Ymaxstd[yi] = Ymaxstd[yi] / Ymax_[yi];         
    }
    
    result.importance = Ymaxstd;
    
    // Find candidate variables that should be added
    // Criteria: a) in principalVars, b) not in currentVars, c) above tolerance
    
    // Create importance-sorted list of speciesPrincipalVars indices only
    result.sortedIndices.setSize(speciesPrincipalVars_.size());
    labelList& sortedIndices = result.sortedIndices;  // Reference for cleaner code
    forAll(sortedIndices, i) { sortedIndices[i] = i; }
    std::sort(sortedIndices.begin(), sortedIndices.end(),
              [&](label a, label b) { 
                  label speciesIdxA = mixture_.species()[speciesPrincipalVars_[a]];
                  label speciesIdxB = mixture_.species()[speciesPrincipalVars_[b]];
                  return Ymaxstd[speciesIdxA] > Ymaxstd[speciesIdxB]; 
              });
    
    Info << nl << "Error-control Variable States:" << nl;
    Info << setw(13) << "Species"
         << setw(18) << "Importance (%)"
         << setw(23) << "Status" << endl;
    Info << "====================================" << nl;

    label numOfNonSpeciesECVars = principalVars_.size() - speciesPrincipalVars_.size();
    result.sortedECIndices.resize(ecVars_.size()-numOfNonSpeciesECVars);
    // Examine each species principalVar in importance order with immediate output
    forAll(sortedIndices, i)
    {
        label speciesIndex = sortedIndices[i];
        word speciesName = speciesPrincipalVars_[speciesIndex];
        // All variables in speciesPrincipalVars_ are guaranteed to be valid species
        label speciesIdx = mixture_.species()[speciesName];
        scalar importance = Ymaxstd[speciesIdx];
        
        // Check criteria using fast lookup operations
        bool notInCurrentVars = !currentVars.found(speciesName);
        bool aboveTolerance = (importance > tolerance);
        
        // Output data immediately during examination
        if (notInCurrentVars && aboveTolerance)
        {
            result.candidateIndices.append(speciesIndex);
            Info << setw(13) << speciesName <<
            setw(18) << fixed << setprecision(2) << importance * 100 <<
            setw(23) << "Candidate" << endl; 
        }
        else if (currentVars.found(speciesName))
        {
            Info << setw(13) << speciesName <<
            setw(18) << fixed << setprecision(2) << importance * 100 <<
            setw(23) << "Already included" << endl;
            result.sortedECIndices.append(speciesIndex);
        }
        else
        {
            Info << setw(13) << speciesName << 
            setw(18) << fixed << setprecision(2) << importance * 100 <<
            setw(23) << "Within tolerance" << endl;
        }
    }
    
    Info << "====================================" << nl << endl;

    // Set underControl flag if optimization is necessary (candidates exist)
    result.underControl = (result.candidateIndices.size() == 0);
    return result;
}


template<class ThermoType>
Foam::scalarField Foam::CCMchemistryModel<ThermoType>::getRRGivenYTP
(
    scalarField& Y,
    scalar& T,
    scalar& p,
    const scalar& deltaT,
    scalar& deltaTChem,
    const scalar& rho,
    const scalar& rho0
)
{
    scalarField Y0(Y);
    scalarField Yupdate(Y);
    scalar Tupdate(T);
    scalar pupdate(p);

    scalar timeLeft = deltaT;
    int dummyCelli=0;

    // Calculate the chemical source terms
    while (timeLeft > small)
    {
        scalar dt = timeLeft;
        solve(pupdate, Tupdate, Yupdate, dummyCelli, dt, deltaTChem);
        timeLeft -= dt;
    }
    
    // update RR
    scalarField RR(nSpecie_);
    for (label i=0; i<nSpecie_; i++)
    {
        RR[i] = (Yupdate[i]*rho - Y0[i]*rho0)/deltaT;
    }
    
    return RR;
}


template<class ThermoType>
Foam::label Foam::CCMchemistryModel<ThermoType>::getNextStepIndex() const
{
    label index = currentStepIndex_;
    currentStepIndex_++;
    
    // Ensure stepTimes_ array is large enough
    if (currentStepIndex_ >= stepTimes_.size())
    {
        // Double the size to accommodate more steps
        stepTimes_.setSize(stepTimes_.size()*2, 0.0);
    }
    
    return index;
}


template<class ThermoType>
typename Foam::CCMchemistryModel<ThermoType>::VariableRanking 
Foam::CCMchemistryModel<ThermoType>::distributeReactionEntry(
    const hashedWordList& encodingVars,
    bool returnRanking
)
{
    // Determine encoding variables
    const hashedWordList& varsToUse = encodingVars.empty() ? principalVars_ : encodingVars;
    
    // Ranges for encoding (either reuse existing or calculate new)
    scalarList encodingMin, encodingMax, encodingSpan;
    // Step 1: Initialization & Setup
    CCM_TIMING_START(Initialization, debugTime_, stepTimer_);
    combustionHelpers_->updatePV();                    // Always update principal variables
    ccmUtilities_->examineRegularVarsMinMax();
    ccmUtilities_->updateMinMaxSpan(varsToUse, encodingMin, encodingMax, encodingSpan, nSlice_);          // Always update principalVars ranges

    tabulation_.reset();

    optionalCpuLoad& chemistryCpuTime
    (
        optionalCpuLoad::New(this->mesh(), "chemistryCpuTime", loadBalancing_)
    );

    // CPU time analysis
    cpuTime solveCpuTime_;

    basicChemistryModel::correct();

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();
    tmp<volScalarField> trho0vf(this->thermo().rho0());
    const volScalarField& rho0vf = trho0vf();

    const volScalarField& T0vf = this->thermo().T().oldTime();
    const volScalarField& p0vf = this->thermo().p().oldTime();

    reactionEvaluationScope scope(*this);

    scalarField Y0(nSpecie_);


    chemistryCpuTime.reset();
    


    // // Initialize reaction entries including Y/T/p and rho
    // List<string>        zoneIndex(mesh().nCells());
    // // labelList       zoneRemainder(mesh().nCells());
    // //----------add for new CCM-------------//
    // labelList       zoneRemainder(mesh().nCells());
    labelList       contactList(RoundRobinCommunicator_.opponents_);
    label contactNum = (Pstream::parRun()) ? RoundRobinCommunicator_.nCores_ : 1;
    
    // Use pre-allocated hash tables and clear them before use
    reactionEntriesList& fromCore = preAllocatedFromCore_;
    reactionEntriesList& toCore = preAllocatedToCore_;
    
    // Clear all hash tables before reuse (19% overhead for 80% performance gain)
    forAll(fromCore, i)
    {
        fromCore[i].clear();
        toCore[i].clear();
    }
    //----------end add for new CCM--------//
    CCM_TIMING_END(Initialization, debugTime_, stepTimer_, *this);

    // logP encoding setup: calculate log(P) min/max/span if enabled
    scalar pMin = -1, pMax = -1, pSpan = -1;
    bool useLogPOnThisRun = true;
    if (highMach_)
    {
        // Calculate log(P) range over all cells
        pMax = gMax(p0vf);
        pMin = max(gMin(p0vf),1);// 1 Pa minimum to avoid log(0)

        Info << nl << "High Mach Pressure: [" << pMin << ", " << pMax << "] Pa, ratio=" << pMax/pMin;

        if (pMax/pMin < 1.1)
        {
            useLogPOnThisRun = false;
            Info << "pMax = " << pMax << ", pMin = " << pMin << ", ratio < 1.1, log(P) encoding disabled" << nl;
        }
        else
        {
            Info << "log(P) encoding enabled" << nl;
            pMax = Foam::log(pMax);
            pMin = Foam::log(pMin);
            pSpan = (pMax - pMin) / nSlice_;
        }
    }

    // Step 2: Cell Grouping & Zone Calculation
    CCM_TIMING_START(Cell_Grouping, debugTime_, stepTimer_);
    groupingTime_ = 0.0;

    forAll(rho0vf, celli)
    {
        const scalar rho = rhovf[celli];
        const scalar rho0 = rho0vf[celli];

        scalar p = p0vf[celli];
        scalar T = T0vf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            Y_[i] = Y0[i] = Yvf_[i].oldTime()[celli];
        }

        // calculate zoneindex and remainder
        string cellZoneIndex = "";


        forAll(varsToUse, vi)
        {
            const volScalarField& field = mesh().lookupObjectRef<volScalarField>(varsToUse[vi]);
            const scalar& value = field[celli];
            label pos = min(max(floor((value-encodingMin[vi])/encodingSpan[vi]),0),maxiRepresentation_);
            cellZoneIndex += parallelComm_().toBase256Word(pos);
        }

        // Add logarithmic pressure suffix if enabled
        string encodedP = "";
        if (useLogPOnThisRun)
        {
            scalar logP = Foam::log(max(p, 1.0));  // Avoid log(0) by clamping to 1 Pa minimum
            label posP = min(max(floor((logP - pMin) / pSpan), 0), maxiRepresentation_);
            string encodedP = parallelComm_().toBase256Word(posP);
            cellZoneIndex += encodedP;
        }
        cellZoneIndex += encodedP;

        //= replace this part
        zoneIndex_[celli] = cellZoneIndex;
        //label remainder = blockWiseRemainder(cellZoneIndex, contactNum);
        uint encodeZoneIndex = parallelComm_().encode(cellZoneIndex);
        label remainder = encodeZoneIndex % contactNum;
        zoneRemainder_[celli] = remainder;


        // Force debug mode if ranking requested
        bool useDebug = debugMode_ || returnRanking;
        ReactionEntry re(Y0, T, p, deltaTChem_[celli], rho0, rho, 1, useDebug);

        groupingTimer_.cpuTimeIncrement();
        reactionEntries& toTable = toCore[remainder];
        if (toTable.found(cellZoneIndex))
        {
            toTable.set(cellZoneIndex, toTable.find(cellZoneIndex)() + re);
        }
        else
        {
            toTable.set(cellZoneIndex, re);
        }
        groupingTime_ += groupingTimer_.cpuTimeIncrement();

    }
    reduce(groupingTime_, maxOp<scalar>());
    CCM_TIMING_END(Cell_Grouping, debugTime_, stepTimer_, *this);

    // Step 3: Round-Robin Distribution
    CCM_TIMING_START(Distribution, debugTime_, stepTimer_);

    // Distribute reactionEntries using round-robin communication
    parallelComm_().performRoundRobinDistribute(toCore, fromCore);

    CCM_TIMING_END(Distribution, debugTime_, stepTimer_, *this);

    CCM_TIMING_START(Merging, debugTime_, stepTimer_);

    gatheredReactionEntries_.clear();
    mergeReactionEntries(gatheredReactionEntries_, fromCore[0]);
    if (Pstream::parRun())
    {
        for (int i=1; i!=contactNum;i++)
        {
            mergeReactionEntries(gatheredReactionEntries_, fromCore[i]);
        }
    }

    CCM_TIMING_END(Merging, debugTime_, stepTimer_, *this);

    // Step 3.3: Post-merge operations (statistics and debug)
    CCM_TIMING_START(Post_merge_Total, debugTime_, stepTimer_);

    label localPhaseSpaceSize = gatheredReactionEntries_.size();
    label phaseSpaceSize      = localPhaseSpaceSize;
    
    CCM_TIMING_END(Post_merge_Total, debugTime_, stepTimer_, *this);
    CCM_TIMING_START(Reduce_PhaseSpace, debugTime_, stepTimer_);
    reduce(phaseSpaceSize, sumOp<label>());
    CCM_TIMING_END(Reduce_PhaseSpace, debugTime_, stepTimer_, *this);
    
    label meanPhaseSpaceSize = phaseSpaceSize/contactNum;
    scalar unbalanceRatio = mag(scalar(localPhaseSpaceSize)/(scalar(meanPhaseSpaceSize)+1)-1.0);
    
    CCM_TIMING_START(Reduce_UnbalanceRatio, debugTime_, stepTimer_);
    reduce(unbalanceRatio, maxOp<scalar>());
    CCM_TIMING_END(Reduce_UnbalanceRatio, debugTime_, stepTimer_, *this);

    // Info << "phaseSpaceSize = " << phaseSpaceSize << endl;
    // Info << "Unbalance ratio = " << unbalanceRatio << endl;
    scalar totalCells = this->mesh().nCells();
    
    CCM_TIMING_START(Reduce_TotalCells, debugTime_, stepTimer_);
    reduce(totalCells, sumOp<scalar>());
    CCM_TIMING_END(Reduce_TotalCells, debugTime_, stepTimer_, *this);

    scalar accRatio = totalCells/scalar(phaseSpaceSize);

    Info << nl << "CCM status" << endl;   
    Info << "====================================" << nl;
    Info << setw(20) << "acceleration ratio"
     << setw(20) << "phase space size"
     << setw(20) << "Unbalance ratio" << endl;

    Info << setw(20) << fixed << setprecision(2) << accRatio
    << setw(20) << fixed << setprecision(2) << phaseSpaceSize
    << setw(20) << fixed << setprecision(2) << unbalanceRatio << endl;
    Info << "====================================" << nl;

    //Info << "acceleration ratio = " << accRatio << endl;
    
    // CCM_TIMING_START(File_Output, debugTime_, stepTimer_);
    // accelerationInfo_.precision(5);
    // accelerationInfo_ << setw(8) << this->mesh().time().value() << tab
    //      << setw(8) << totalCells << tab
    //      << setw(8) << phaseSpaceSize << tab
    //      << setw(8) << accRatio << endl;
    // CCM_TIMING_END(File_Output, debugTime_, stepTimer_, *this);

    CCM_TIMING_START(Debug_Operations, debugTime_, stepTimer_);
    if (debugMode_)
    {
        ccmUtilities_->updateYMax();
        ccmDebug_->debugStdRank(gatheredReactionEntries_);
        ccmDebug_->debugAcc4NonTrivalCells(gatheredReactionEntries_);
    }
    CCM_TIMING_END(Debug_Operations, debugTime_, stepTimer_, *this);

    // Return ranking if requested
    VariableRanking result;
    if (returnRanking)
    {
        ccmUtilities_->updateYMax();  // Ensure Ymax is current
        result = this->calculateVariableRanking(encodingVars.empty() ? principalVars_ : encodingVars, 1.0/nSlice_);
        
        if (debugMode_)
        {
            ccmDebug_->debugStdRank(gatheredReactionEntries_);  // Still output debug info
        }
    }
    
    return result;
}

template<class ThermoType>
Foam::hashedWordList 
Foam::CCMchemistryModel<ThermoType>::findLeastImportantEcVars(
    const VariableRanking& ranking, 
    label n
) const
{
    hashedWordList varsToRemove;
    
    // Safety check: ensure we don't remove all variables (keep at least 1)
    label maxToRemove = min(n, max(0, label(ranking.sortedECIndices.size() - 1)));
    
    if (maxToRemove <= 0)
    {
        Info << "Warning: Cannot remove variables - ecVars has " 
             << ranking.sortedECIndices.size() << " variables, keeping at least 1" << endl;
        return varsToRemove;
    }
    
    // Take the last n elements (least important) from sortedECIndices
    // Since sortedECIndices is already sorted by importance (descending)
    for (label i = ranking.sortedECIndices.size() - maxToRemove; 
         i < ranking.sortedECIndices.size(); i++)
    {
        word varName = speciesPrincipalVars_[ranking.sortedECIndices[i]];
        varsToRemove.append(varName);
    }
    
    return varsToRemove;
}

template<class ThermoType>
void Foam::CCMchemistryModel<ThermoType>::updateReactionRate()
{

    // Step 4: Chemistry Computation
    CCM_TIMING_START(Chemistry_Computation, debugTime_, stepTimer_);
    // Calculate reaction rates for the reactionEntries assigned to this processor 
    // according to the remainders
    // Store the results in the reaction rate entries(rrEntries)


    // Use pre-allocated local RR entries and clear before use
    rrEntries& localRREntries = preAllocatedLocalRREntries_;
    localRREntries.clear();

    for (auto it = gatheredReactionEntries_.begin(); it != gatheredReactionEntries_.end(); it++)
    {
        ReactionEntry re = it();
        scalarField& Y = re.Y;
        scalar& T = re.T;
        scalar& p = re.p;
        const scalar& deltaT = mesh().time().deltaTValue();
        scalar deltaTChem = re.dtChem;
        const scalar& rho0 = re.rho0;
        const scalar& rho = re.rho;

        scalarField cellRR = getRRGivenYTP(Y, T, p, deltaT, deltaTChem, rho, rho0);
        if (debugMode_)
        {
            localRREntries.set(it.key(), RREntry(cellRR, deltaTChem,debugMode_,re.Y));
        }
        else
        {
            localRREntries.set(it.key(), RREntry(cellRR, deltaTChem));
        }
    }

    CCM_TIMING_END(Chemistry_Computation, debugTime_, stepTimer_, *this);
}

template<class ThermoType>
void Foam::CCMchemistryModel<ThermoType>::distributeReactionRate()
{

    reactionEntriesList& fromCore = preAllocatedFromCore_;
    label contactNum = (Pstream::parRun()) ? RoundRobinCommunicator_.nCores_ : 1;
    rrEntries& localRREntries = preAllocatedLocalRREntries_;
    rrEntriesList& globalRREntries(preAllocatedReceivedRREntries_);
    // Step 5: Communication Infrastructure Setup
    CCM_TIMING_START(Comm_Setup, debugTime_, stepTimer_);
    // Choose communication method based on configuration
    //rrEntriesList globalRREntries(RoundRobinCommunicator_.nCores_);
    CCM_TIMING_END(Comm_Setup, debugTime_, stepTimer_, *this);

    if (optimizedCommunication_)
    {
        // Step 6: Targeted Hash Table Construction (Optimized method only)
        CCM_TIMING_START(Hash_Table_Construction, debugTime_, stepTimer_);
        // NEW OPTIMIZED METHOD: Targeted distribution using pre-allocated hash tables
        rrEntriesList& returnToCore = preAllocatedReturnToCore_;
        
        // Clear pre-allocated communication hash tables
        forAll(returnToCore, i)
        {
            returnToCore[i].clear();
        }
        
        // Create targeted return data by looking up results
        for (int coreIndex = 0; coreIndex < contactNum; coreIndex++)
        {
            reactionEntries& zonesFromThisCore = fromCore[coreIndex];
            rrEntries& resultsForThisCore = returnToCore[coreIndex];
            
            for (auto it = zonesFromThisCore.begin(); it != zonesFromThisCore.end(); it++)
            {
                string zoneIndex = it.key();
                if (localRREntries.found(zoneIndex))
                {
                    resultsForThisCore.set(zoneIndex, localRREntries.find(zoneIndex)());
                }
            }
        }
        CCM_TIMING_END(Hash_Table_Construction, debugTime_, stepTimer_, *this);

        // Step 7: Communication Distribution/Broadcast
        CCM_TIMING_START(Communication_Execution, debugTime_, stepTimer_);
        // Targeted distribution using pre-allocated hash tables
        rrEntriesList& receivedRREntries = preAllocatedReceivedRREntries_;
        
        // Clear received RR entries before use
        // forAll(receivedRREntries, i)
        // {
        //     receivedRREntries[i].clear();
        // }
        
        // Use parallelComm_().performRoundRobinDistribute function
        parallelComm_().performRoundRobinDistribute(returnToCore, receivedRREntries);
        
    }
    else
    {
        // Step 6: No hash table construction for original method
        // Step 6 time is automatically 0 for original method - no manual adjustment needed
        // Step 7: Communication Distribution/Broadcast
        CCM_TIMING_START(Communication_Execution, debugTime_, stepTimer_);
        
        // ORIGINAL METHOD: Broadcast all results to all cores
        label myCoreNumLCS = RoundRobinCommunicator_.myCoreNumLCS_;
        parallelComm_().performRoundRobinBroadcast(localRREntries, globalRREntries, myCoreNumLCS);
    }
    CCM_TIMING_END(Communication_Execution, debugTime_, stepTimer_, *this);

    // // Step 8: Rate Assignment to Cells
    //rrEntriesList& globalRREntries = preAllocatedReceivedRREntries_;
    CCM_TIMING_START(Rate_Assignment, debugTime_, stepTimer_);
    // // Attribute the reaction rates and the deltaChemT
    forAll(zoneRemainder_, celli)
    {
        const label& remainder = zoneRemainder_[celli];
        rrEntries& lookupFrom = globalRREntries[remainder];
        scalarField cellRR = lookupFrom.find(zoneIndex_[celli])().RR;
        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = cellRR[i];
        }
        deltaTChem_[celli] = lookupFrom.find(zoneIndex_[celli])().deltaTChem;
        deltaTChem_[celli] = min(deltaTChem_[celli], deltaTChemMax_);

        if (debugMode_ && examineYdiff_)
        {
            const scalarField& Yavg = lookupFrom.find(zoneIndex_[celli])().Y;
            for (label i=0; i<nSpecie_; i++)
            {
                Ydiff_[i][celli] = Yavg[i] - Yvf_[i].oldTime()[celli];
                Ydiff_[i][celli] /= Ymax_[i];
            }
        }
    }
    if (debugMode_ && examineYdiff_ && shutdownImmediately_)
    {
        // write out all Ydiff_
        forAll(Ydiff_, i)
        {
            Ydiff_[i].write();
        }
        SHUTDOWN
    }
    //deltaTMin = min(deltaTChem_).value();
    CCM_TIMING_END(Rate_Assignment, debugTime_, stepTimer_, *this);
}

template<class ThermoType>
void Foam::CCMchemistryModel<ThermoType>::updateCCM4MeshChange()
{
    if (mesh().dynamic())
    {
        label nCells = mesh().nCells();
        if (zoneIndex_.size() != nCells)
        {
           
            zoneIndex_.setSize(nCells);
            zoneRemainder_.setSize(nCells);
        }
    }
}






// ************************************************************************* //
