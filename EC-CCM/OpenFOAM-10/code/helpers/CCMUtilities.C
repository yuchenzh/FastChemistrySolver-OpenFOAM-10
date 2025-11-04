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

#include "CCMUtilities.H"
#include "atomicWeights.H"
#include "ops.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::CCMUtilities<ThermoType>::CCMUtilities
(
    const multiComponentMixture<ThermoType>& mixture,
    const PtrList<ThermoType>& specieThermos,
    const dictionary& CCMdict,
    const label& nSpecie,
    hashedWordList& principalVars,
    hashedWordList& updateVars,
    hashedWordList& regularVars,
    const Switch& oldStylePhi,
    scalarList& JhCoeff,
    scalarList& JcCoeff,
    scalarList& JnCoeff,
    scalarList& JoCoeff,
    scalarList& Jh_h2oCoeff,
    scalarList& Jc_co2Coeff,
    scalarList& Jo_co2_h2oCoeff,
    const PtrList<volScalarField>& Yvf,
    scalarField& Ymax,
    HashTable<scalar, word>& pvMin,
    HashTable<scalar, word>& pvMax,
    HashTable<scalar, word>& pvSpan,
    const scalar& ignoreMin,
    const fvMesh& mesh
)
:
    mixture_(mixture),
    specieThermos_(specieThermos),
    CCMdict_(CCMdict),
    nSpecie_(nSpecie),
    principalVars_(principalVars),
    updateVars_(updateVars),
    regularVars_(regularVars),
    oldStylePhi_(oldStylePhi),
    JhCoeff_(JhCoeff),
    JcCoeff_(JcCoeff),
    JnCoeff_(JnCoeff),
    JoCoeff_(JoCoeff),
    Jh_h2oCoeff_(Jh_h2oCoeff),
    Jc_co2Coeff_(Jc_co2Coeff),
    Jo_co2_h2oCoeff_(Jo_co2_h2oCoeff),
    Yvf_(Yvf),
    Ymax_(Ymax),
    pvMin_(pvMin),
    pvMax_(pvMax),
    pvSpan_(pvSpan),
    ignoreMin_(ignoreMin),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::CCMUtilities<ThermoType>::appendIfNotExisting
(
    const Foam::word& varName, 
    Foam::hashedWordList& listToBeAppended
)
{
    bool found = false;
    found = listToBeAppended.found(varName);
    if (!found)
    {
        listToBeAppended.append(varName);
    }
    return found;
}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::outputHashedWordList(const hashedWordList& list) const
{
    forAll(list, i)
    {
        Info << list[i] << " ";
    }
    Info << endl;
}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::initPvLists()
{
    // Deal with variables to update
    if (principalVars_.found("phieq"))
    {
        if (oldStylePhi_)
        {
            appendIfNotExisting("Jc", updateVars_);
            appendIfNotExisting("Jh", updateVars_);
            appendIfNotExisting("Jo", updateVars_);
            appendIfNotExisting("phieq", updateVars_);
        }
        else
        {
            appendIfNotExisting("phieq", updateVars_);
            appendIfNotExisting("Jc_co2", updateVars_);
            appendIfNotExisting("Jh_h2o", updateVars_);
            appendIfNotExisting("Jo_co2_h2o", updateVars_);
        }
    }
    if (principalVars_.found("chi"))
    {
        appendIfNotExisting("chi", updateVars_);
        if (!CCMdict_.found("chiComponent"))
        {
            FatalErrorIn("CCMUtilities<ThermoType>::initPvLists()")
                << "The chiComponent entry not found in CCM dictionary" 
                << exit(FatalError);
        }
    }
    if (principalVars_.found("J"))
    {
        appendIfNotExisting("J", updateVars_);
        if (!CCMdict_.found("JElement"))
        {
            FatalErrorIn("CCMUtilities<ThermoType>::initPvLists()")
                << "The JElement entry not found in CCM dictionary"
                << exit(FatalError);
        }
        else
        {
            word JElement = CCMdict_.lookup("JElement");
            if (JElement == "C" || JElement == "c")
            {
                appendIfNotExisting("Jc", updateVars_);
            }
            else if (JElement == "H" || JElement == "h")
            {
                appendIfNotExisting("Jh", updateVars_);
            }
            else if (JElement == "N" || JElement == "n")
            {
                appendIfNotExisting("Jn", updateVars_);
            }
            else if (JElement == "O" || JElement == "o")
            {
                appendIfNotExisting("Jo", updateVars_);
            }
            else
            {
                FatalErrorIn("CCMUtilities<ThermoType>::initPvLists()")
                    << "The JElement entry should be one of C, H, N, O"
                    << exit(FatalError);
            }
        }
    }

    // Deal with species-only principal variables
    hashedWordList speciesPrincipalVars;
    forAll(principalVars_, pvi)
    {
        forAll(mixture_.species(), i)
        {
            if (mixture_.species()[i] == principalVars_[pvi])
            {
                speciesPrincipalVars.append(principalVars_[pvi]);
                break;
            }
        }
    }

    Info << "Principal variables: ";
    outputHashedWordList(principalVars_);
    Info << "Regular variables: ";
    outputHashedWordList(regularVars_);
    Info << "Variables to be update: ";
    outputHashedWordList(updateVars_);
}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::initRegularVars()
{
    dictionary pvDict(CCMdict_.subDict("pvInfo"));
    List<Foam::keyType> keys = pvDict.keys();
    regularVars_.clear();

    forAll(keys, ki)
    {
        appendIfNotExisting(keys[ki], regularVars_);
    }
}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::initJCoeffs()
{
    label H2Oindex(-1);
    label CO2index(-1);

    if (mixture_.species().found("H2O"))
    {
        H2Oindex = mixture_.species()["H2O"];
    }
    if (mixture_.species().found("CO2"))
    {
        CO2index = mixture_.species()["CO2"];
    }
    
    // Initialize JcCoeff_, JnCoeff_, JoCoeff_, JhCoeff_
    const scalar MWh = atomicWeights["H"];
    const scalar MWc = atomicWeights["C"];
    const scalar MWn = atomicWeights["N"];
    const scalar MWo = atomicWeights["O"];

    for (label i=0; i<nSpecie_; i++)
    {
        const List<specieElement>& spComposition = mixture_.specieComposition(i);
        forAll(spComposition, elementI)
        {
            scalar MW = specieThermos_[i].W();
            if (spComposition[elementI].name() == "H")
            {
                JhCoeff_[i] = spComposition[elementI].nAtoms() * MWh / MW;
            }
            if (spComposition[elementI].name() == "C")
            {
                JcCoeff_[i] = spComposition[elementI].nAtoms() * MWc / MW;
            }
            if (spComposition[elementI].name() == "N")
            {
                JnCoeff_[i] = spComposition[elementI].nAtoms() * MWn / MW;
            }
            if (spComposition[elementI].name() == "O")
            {
                JoCoeff_[i] = spComposition[elementI].nAtoms() * MWo / MW;
            }
        }
    }

    // Initialize modified J coefficients
    Jh_h2oCoeff_ = JhCoeff_;
    Jc_co2Coeff_ = JcCoeff_;
    Jo_co2_h2oCoeff_ = JoCoeff_;

    if (H2Oindex != -1)
    {
        Jh_h2oCoeff_[H2Oindex] = 0.0;
        Jo_co2_h2oCoeff_[H2Oindex] = 0.0;
    }

    if (CO2index != -1)
    {
        Jc_co2Coeff_[CO2index] = 0.0;
        Jo_co2_h2oCoeff_[CO2index] = 0.0;
    }

}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::examineMinMax() const
{

    Info << "\nCCM Variable Analysis:" << endl;
    Info << "====================================" << endl;
    Info << setw(15) << "Variable" << setw(12) << "Min Value" << setw(12) << "Max Value" << endl;
    
    forAll(principalVars_, pvi)
    {
        const word& varName = principalVars_[pvi];
        const volScalarField& field = mesh_.lookupObjectRef<volScalarField>(varName);
        Info << setw(15) << varName 
            << setw(12) << scientific << setprecision(4) << gMin(field)
            << setw(12) << scientific << setprecision(4) << gMax(field)
            << endl;
    }
    Info << "#====== CCM Information ======#" << endl;

}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::examineRegularVarsMinMax() const
{
    if (regularVars_.size() >0)
    {
        Info << "\nCCM Regular Variable Analysis:" << endl;
        Info << "====================================" << endl;
        Info << setw(15) << "Variable" << setw(12) << "Min Value" << setw(12) << "Max Value" << endl;
        
        forAll(regularVars_, ri)
        {
            const word& varName = regularVars_[ri];
            // loookup the field and get its min/max value
            const volScalarField& field = mesh_.lookupObjectRef<volScalarField>(varName);
            Info << setw(15) << varName 
                << setw(12) << scientific << setprecision(4) << gMin(field)
                << setw(12) << scientific << setprecision(4) << gMax(field)
                << endl;
        }
        Info << "#====== CCM Information ======#" << endl;
    }
}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::updateYMax()
{
    forAll(Yvf_, i)
    {
        scalar tempYmax = gMax(Yvf_[i]);
        scalar threshold = 1e-8;
        Ymax_[i] = (tempYmax > threshold) ? tempYmax : 1.0;
    }
    reduce(Ymax_, maxOp<scalarList>());
}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::initPvMinMaxSpan()
{
    // Check whether regularVars_ is a subset of principalVars_
    forAll (regularVars_, ri)
    {
        if (!principalVars_.found(regularVars_[ri]))
        {
            FatalErrorIn("CCMUtilities<ThermoType>::initPvMinMaxSpan()")
                << "The regular variable " << regularVars_[ri] 
                << " is not a subset of principalVars_" << exit(FatalError);
        }
    }

    // Read pvInfo from dictionary
    dictionary pvDict(CCMdict_.subDict("pvInfo"));
    forAll(regularVars_, ri)
    {
        word varName = regularVars_[ri];
        List<scalar> values = pvDict.lookup(varName);
        if (values.size() != 3)
        {
            FatalErrorIn("CCMUtilities<ThermoType>::initPvMinMaxSpan()")
                << "The pvInfo should be in the format \"varName (min max span);\""
                << exit(FatalError);
        }
        pvMin_.insert(varName, values[0]);
        pvMax_.insert(varName, values[1]);
        pvSpan_.insert(varName, values[2]);
    }

    // Initialize remaining variables with default values
    forAll(principalVars_, pvi)
    {
        word varName = principalVars_[pvi];
        if (!regularVars_.found(varName))
        {
            pvMin_.insert(varName, 0.0);
            pvMax_.insert(varName, 1.0);
            pvSpan_.insert(varName, ignoreMin_/10.0);
        }
    }
}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::updateMinMaxSpan(
    const hashedWordList& vars, 
    scalarList& minVals, 
    scalarList& maxVals, 
    scalarList& spanVals,
    const label& nSlice
) const
{
    minVals.resize(vars.size());
    maxVals.resize(vars.size());
    spanVals.resize(vars.size());
    
    forAll(vars, vi)
    {
        const word& varName = vars[vi];
        if (regularVars_.found(varName))
        {
            // Use configured values from hash tables - direct O(1) access
            minVals[vi] = pvMin_[varName];
            maxVals[vi] = pvMax_[varName];
            spanVals[vi] = pvSpan_[varName];
        }
        else
        {
            // Auto-calculate (same logic as original)
            const volScalarField& field = mesh_.lookupObjectRef<volScalarField>(varName);
            minVals[vi] = gMin(field);//max(gMin(field), 0.);
            maxVals[vi] = max(gMax(field), minVals[vi] + ignoreMin_);
            spanVals[vi] = max((maxVals[vi] - minVals[vi]) / nSlice, ignoreMin_/10.);
        }
    }
}


template<class ThermoType>
void Foam::CCMUtilities<ThermoType>::updatePvMinMaxSpan(const label& nSlice)
{
    scalarList minVals, maxVals, spanVals;
    updateMinMaxSpan(principalVars_, minVals, maxVals, spanVals, nSlice);
    
    forAll(principalVars_, pi)
    {
        const word& pvName = principalVars_[pi];
        if (!regularVars_.found(pvName))
        {
            pvMin_.insert(pvName, minVals[pi]);
            pvMax_.insert(pvName, maxVals[pi]);
            pvSpan_.insert(pvName, spanVals[pi]);
        }
    }
}


// ************************************************************************* //