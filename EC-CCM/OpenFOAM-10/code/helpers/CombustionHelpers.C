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

#include "CombustionHelpers.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::CombustionHelpers<ThermoType>::CombustionHelpers
(
    const fvMesh& mesh,
    const fluidReactionThermo& thermo,
    const dictionary& CCMdict,
    const PtrList<volScalarField>& Yvf,
    const hashedWordList& updateVars,
    const scalarList& JhCoeff,
    const scalarList& JcCoeff,
    const scalarList& JnCoeff,
    const scalarList& JoCoeff,
    const scalarList& Jh_h2oCoeff,
    const scalarList& Jc_co2Coeff,
    const scalarList& Jo_co2_h2oCoeff,
    volScalarField& Jh,
    volScalarField& Jc,
    volScalarField& Jn,
    volScalarField& Jo,
    volScalarField& Jh_h2o,
    volScalarField& Jc_co2,
    volScalarField& Jo_co2_h2o,
    volScalarField& J,
    volScalarField& phieq,
    volScalarField& chi,
    const Switch& oldStyleChi,
    const Switch& oldStylePhi,
    const scalar& fuelOtoC
)
:
    mesh_(mesh),
    thermo_(thermo),
    CCMdict_(CCMdict),
    Yvf_(Yvf),
    updateVars_(updateVars),
    JhCoeff_(JhCoeff),
    JcCoeff_(JcCoeff),
    JnCoeff_(JnCoeff),
    JoCoeff_(JoCoeff),
    Jh_h2oCoeff_(Jh_h2oCoeff),
    Jc_co2Coeff_(Jc_co2Coeff),
    Jo_co2_h2oCoeff_(Jo_co2_h2oCoeff),
    Jh_(Jh),
    Jc_(Jc),
    Jn_(Jn),
    Jo_(Jo),
    Jh_h2o_(Jh_h2o),
    Jc_co2_(Jc_co2),
    Jo_co2_h2o_(Jo_co2_h2o),
    J_(J),
    phieq_(phieq),
    chi_(chi),
    oldStyleChi_(oldStyleChi),
    oldStylePhi_(oldStylePhi),
    fuelOtoC_(fuelOtoC)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJc()
{
    Jc_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jc_ += JcCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJh()
{
    Jh_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jh_ += JhCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJo()
{
    Jo_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jo_ += JoCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJn()
{
    Jn_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jn_ += JnCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJh_h2o()
{
    Jh_h2o_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jh_h2o_ += Jh_h2oCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJc_co2()
{
    Jc_co2_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jc_co2_ += Jc_co2Coeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJo_co2_h2o()
{
    Jo_co2_h2o_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jo_co2_h2o_ += Jo_co2_h2oCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateChi()
{
    if (oldStyleChi_)
    {
        volVectorField gradMZ_0 = fvc::grad(phieq_);
        dimensionedScalar dimDiffutionCoeff(dimLength*dimLength/dimTime, 1.);
        chi_ = dimDiffutionCoeff*(gradMZ_0&gradMZ_0);
        forAll(chi_, i)
        {
            chi_[i] = Foam::exp(-chi_[i]);
        }
        forAll(chi_.boundaryField(), patchi)
        {
            forAll(chi_.boundaryField()[patchi], facei)
            {
                chi_.boundaryFieldRef()[patchi][facei] = Foam::exp(-chi_.boundaryField()[patchi][facei]);
            }
        }
    }
    else
    {
        volScalarField D
        (
            mesh_.lookupObjectRef<volScalarField>("thermo:kappa")
            /
            thermo_.Cp()/thermo_.rho()
        );
        const word& chiComponent = CCMdict_.lookup("chiComponent");
        volScalarField& Z(mesh_.objectRegistry::lookupObjectRef<volScalarField>(chiComponent));
        volVectorField gradZ(fvc::grad(Z));
        chi_ = 2*D*(gradZ&gradZ);
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updatePhi()
{
    if (oldStylePhi_)
    {
        phieq_ = 2*Jc_/atomicWeights["C"] + 0.5*Jh_/atomicWeights["H"]-fuelOtoC_*Jc_/atomicWeights["C"];
        phieq_ /= max(Jo_/atomicWeights["O"]-fuelOtoC_*Jc_/atomicWeights["C"], 1e-4);
    }
    else
    {
        phieq_ = 2*Jc_co2_/atomicWeights["C"] + 0.5*Jh_h2o_/atomicWeights["H"]
        -fuelOtoC_*Jc_co2_/atomicWeights["C"];
        phieq_ /= max(Jo_co2_h2o_/atomicWeights["O"]-fuelOtoC_*Jc_co2_/atomicWeights["C"], 1e-4);        
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updatePV()
{
    // Always update J* first because phieq and chi depend on them
    if (updateVars_.found("Jc"))
    {
        updateJc();
    }
    if (updateVars_.found("Jh"))
    {
        updateJh();
    }
    if (updateVars_.found("Jo"))
    {
        updateJo();
    }
    if (updateVars_.found("Jn"))
    {
        updateJn();
    }
    if (updateVars_.found("Jh_h2o"))
    {
        updateJh_h2o();
    }
    if (updateVars_.found("Jc_co2"))
    {
        updateJc_co2();
    }
    if (updateVars_.found("Jo_co2_h2o"))
    {
        updateJo_co2_h2o();
    }
    
    // Update J based on JElement setting
    if (updateVars_.found("J"))
    {
        const word& JElement = CCMdict_.lookup("JElement");
        if (JElement == "h")
        {
            J_ = Jh_;
        }
        else if (JElement == "c")
        {
            J_ = Jc_;
        }
        else if (JElement == "o")
        {
            J_ = Jo_;
        }
        else if (JElement == "n")
        {
            J_ = Jn_;
        }
        else
        {
            FatalErrorIn("CombustionHelpers<ThermoType>::updatePV()")
                << "Unknown JElement: " << JElement << exit(FatalError);
        }
    }
    
    // Update phieq (depends on J*)
    if (updateVars_.found("phieq"))
    {
        updatePhi();
    }
    
    // Update chi (depends on other variables)
    if (updateVars_.found("chi"))
    {
        updateChi();
    }
}


// ************************************************************************* //