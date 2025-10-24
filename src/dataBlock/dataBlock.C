/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::dataBlock

Description

SourceFiles
    dataBlock.C

\*---------------------------------------------------------------------------*/

#include "dataBlock.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dataBlock::dataBlock(label& nSpecies_):
Y_(nSpecies_),
T(300.0),
p(101325),
CPUtime(0),
celli(-1),
deltaT(0),
deltaTChem_(0)
{}

Foam::dataBlock::dataBlock(Istream& is)
{
    this->Y_.clear();
    is >> this->Y_;
    is >> this->T;
    is >> this->p;
    is >> this->CPUtime;
    is >> this->celli;
    is >> this->deltaT;
    is >> this->deltaTChem_;
}

// * * * * * * * * * * * * * * * IOstream Operator  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator >>(Istream& is, dataBlock& data)
{
    is >> data.Y_;
    is >> data.T;
    is >> data.p;
    is >> data.CPUtime;
    is >> data.celli;
    is >> data.deltaT;
    is >> data.deltaTChem_;
    return is;
}

Foam::Ostream& Foam::operator <<(Ostream& os, const dataBlock& data)
{
    os << data.Y_;
    os << data.T;
    os << data.p;
    os << data.CPUtime;
    os << data.celli;
    os << data.deltaT;
    os << data.deltaTChem_;
    return os;
}






// ************************************************************************* //
