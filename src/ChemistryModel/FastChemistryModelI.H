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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline label
Foam::FastChemistryModel<ThermoType>::nSpecie() const
{
    return nSpecie_;
}

template<class ThermoType>
inline Foam::PtrList<Foam::volScalarField>&
Foam::FastChemistryModel<ThermoType>::RR()
{
    return RR_;
}

template<class ThermoType>
inline Foam::label Foam::FastChemistryModel<ThermoType>::nReaction() const
{
    return reaction.n_Reactions;
}

template<class ThermoType>
inline OptReaction& Foam::FastChemistryModel<ThermoType>::getReaction() const
{
    return reaction;
}


template<class ThermoType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FastChemistryModel<ThermoType>::RR(const label i) const
{
    return RR_[i].internalField();
}


template<class ThermoType>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FastChemistryModel<ThermoType>::RR(const label i)
{
    return RR_[i].ref();
}




// ************************************************************************* //
