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

#include "CCMDebug.H"
#include "ops.H"
#include "Pstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::CCMDebug<ThermoType>::CCMDebug
(
    const label& nSpecie,
    const scalarField& Ymax,
    const multiComponentMixture<ThermoType>& mixture,
    scalarField& stepTimes,
    wordList& stepNames,
    label& currentStepIndex
)
:
    nSpecie_(nSpecie),
    Ymax_(Ymax),
    mixture_(mixture),
    stepTimes_(stepTimes),
    stepNames_(stepNames),
    currentStepIndex_(currentStepIndex)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::CCMDebug<ThermoType>::debugStdRank
(
    const reactionEntries& rEs
) const
{
    scalarField Ymaxstd(nSpecie_, 0.0);
    label i = 0;
    for (auto it = rEs.begin(); it != rEs.end(); it++)
    {
        const ReactionEntry& re = it();
        forAll(re.Y, yi)
        {
            if (re.count > 1)
            {
                Ymaxstd[yi] = max(Ymaxstd[yi], re.Ystd[yi]);
            }
        }
        i++;
    }
    i = 0;
    // deal with Ystd

    reduce(Ymaxstd, maxOp<scalarList>());
    // calculate relevant importance
    forAll(Ymax_, yi)
    {
        Ymaxstd[yi] = Ymaxstd[yi] / (Ymax_[yi]);         
    }

    // now let's sort the Ymaxstd, and output the rank
    List<std::pair<scalar, label> > YmaxstdPairs(nSpecie_);
    forAll(Ymaxstd, yi)
    {
        YmaxstdPairs[yi] = std::make_pair(Ymaxstd[yi], yi);
    }
    
    // Sort in descending order
    std::sort(YmaxstdPairs.begin(), YmaxstdPairs.end(), 
              [](const std::pair<scalar, label>& a, const std::pair<scalar, label>& b) 
              {
                  return a.first > b.first;
              });

    Info << "Species ranking by normalized standard deviation (only species with norm_std > 1e-10 are displayed):" << endl;
    Info << "Rank\tSpecies\t\tNorm_Std" << endl;
    forAll(YmaxstdPairs, i)
    {
        if (YmaxstdPairs[i].first > 1e-10)  // Only show significant entries
        {
            label speciesIndex = YmaxstdPairs[i].second;
            Info << i+1 << "\t" << mixture_.species()[speciesIndex] 
                 << "\t\t" << YmaxstdPairs[i].first << endl;
        }
    }
}


template<class ThermoType>
void Foam::CCMDebug<ThermoType>::debugAcc4NonTrivalCells
(
    const reactionEntries& rEs
) const
{
    label localPhaseSpaceSize = rEs.size();
    
    label localNonTrivialCells = 0;
    label localTotalCells = 0;
    
    for (auto it = rEs.begin(); it != rEs.end(); it++)
    {
        const ReactionEntry& re = it();
        localTotalCells += re.count;
        if (re.count > 1)
        {
            localNonTrivialCells += re.count;
        }
    }
    
    // Reduce across all processors
    label globalPhaseSpaceSize = localPhaseSpaceSize;
    label globalNonTrivialCells = localNonTrivialCells;
    label globalTotalCells = localTotalCells;
    
    reduce(globalPhaseSpaceSize, sumOp<label>());
    reduce(globalNonTrivialCells, sumOp<label>());
    reduce(globalTotalCells, sumOp<label>());
    
    scalar accelerationRatio = 0.0;
    if (globalPhaseSpaceSize > 0)
    {
        accelerationRatio = scalar(globalTotalCells) / scalar(globalPhaseSpaceSize);
    }
    
    if (Pstream::master())
    {
        Info << "phaseSpaceSize = " << globalPhaseSpaceSize << endl;
        Info << "Unbalance ratio = " << scientific << setprecision(4) 
             << scalar(globalNonTrivialCells) / scalar(globalTotalCells) << endl;
        Info << "acceleration ratio = " << scientific << setprecision(4) 
             << accelerationRatio << endl;
    }
}


template<class ThermoType>
void Foam::CCMDebug<ThermoType>::outputTimingAnalysis
(
    const scalarField& stepTimes
) const
{
    // Calculate total time
    scalar totalTime = sum(stepTimes);
    
    // Calculate percentages based on actual executed steps
    label actualSteps = stepNames_.size();
    scalarField percentages(actualSteps);
    for (label i = 0; i < actualSteps; i++)
    {
        percentages[i] = (totalTime > 0) ? 100.0 * stepTimes[i] / totalTime : 0.0;
    }

    // Reduce timing data across all processors for parallel runs - only for actual used steps
    scalarField globalStepTimes(actualSteps, 0.0);
    scalarField maxStepTimes(actualSteps, 0.0);
    scalarField minStepTimes(actualSteps, 0.0);
    
    // Copy actual data
    for (label i = 0; i < actualSteps; i++)
    {
        globalStepTimes[i] = stepTimes[i];
        maxStepTimes[i] = stepTimes[i];
        minStepTimes[i] = stepTimes[i];
    }
    
    if (Pstream::parRun())
    {
        reduce(globalStepTimes, sumOp<scalarField>());
        reduce(maxStepTimes, maxOp<scalarField>());
        reduce(minStepTimes, minOp<scalarField>());
        
        // Calculate average times
        globalStepTimes /= Pstream::nProcs();
    }

    // Output detailed timing information (only on master processor)
    if (Pstream::master())
    {
        Info << nl << "CCM Solve Function Timing Analysis:" << nl;
        Info << "====================================" << nl;
        
        // Header
        Info << setw(25) << "Step" << setw(12) << "Time(s)" 
             << setw(10) << "Percent";
        
        if (Pstream::parRun())
        {
            Info << setw(12) << "Max(s)" << setw(12) << "Min(s)";
        }
        Info << nl;
        
        // Data rows
        for (label i = 0; i < actualSteps; i++)
        {
            Info << setw(25) << stepNames_[i]
                 << setw(12) << fixed << setprecision(6) << globalStepTimes[i]
                 << setw(9) << fixed << setprecision(2) << percentages[i] << "%";
            
            if (Pstream::parRun())
            {
                Info << setw(12) << fixed << setprecision(6) << maxStepTimes[i]
                     << setw(12) << fixed << setprecision(6) << minStepTimes[i];
            }
            Info << nl;
        }
        
        Info << "====================================" << nl;
        Info << setw(25) << "Total"
             << setw(12) << fixed << setprecision(6) << totalTime
             << setw(9) << "100.00%" << nl << endl;;
    }
}


template<class ThermoType>
Foam::label Foam::CCMDebug<ThermoType>::getNextStepIndex() const
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


// ************************************************************************* //