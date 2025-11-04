/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#ifndef REACTION_ENTRY_H
#define REACTION_ENTRY_H

#include "fvCFD.H"

namespace Foam
{

// declaration
Foam::word num2Word(const Foam::label& num, const Foam::label& digits);

class ReactionEntry;
Ostream& operator << (Ostream& os, const ReactionEntry& re);

Istream& operator >> (Istream& is, ReactionEntry& re);

bool operator == (const ReactionEntry& re1, const ReactionEntry& re2);

bool operator != (const ReactionEntry& re1, const ReactionEntry& re2);

ReactionEntry operator+ (const ReactionEntry& re1, const ReactionEntry& re2);

using reactionEntries = HashTable<ReactionEntry, string>;
using reactionEntriesList = List<reactionEntries>;
using reactionEntriesLists = List<reactionEntriesList>;

void mergeReactionEntries(reactionEntries& re1, const reactionEntries& re2);


class ReactionEntry
{

    public:
    //- constructor
    ReactionEntry
    (
        const scalarField& Yvalue,
        const scalar& Tvalue,
        const scalar& pValue,
        const scalar& dtChemValue,
        const scalar& rho0Value,
        const scalar& rhoValue,
        const label& countValue,
        bool debugValue = false,
        const scalarField& YstdValue = scalarField(0)
    ):
    Y(Yvalue), 
    T(Tvalue), 
    p(pValue), 
    dtChem(dtChemValue), 
    rho0(rho0Value),
    rho(rhoValue),
    count(countValue),
    Ystd(YstdValue),
    debug(debugValue)
    {}

    // // construct from Foam::Istream&
    ReactionEntry(Istream& is)
    {
        read(is);
    }

    // //  default constructor
    ReactionEntry():
    Y(scalarField(0)),
    T(0),
    p(0),
    dtChem(0),
    rho0(0),
    rho(0),
    count(0),
    Ystd(scalarField(0)),
    debug(false)
    {}

    // // copy constructor
    ReactionEntry(const ReactionEntry& re):
    Y(re.Y),
    T(re.T),
    p(re.p),
    dtChem(re.dtChem),
    rho0(re.rho0),
    rho(re.rho),
    count(re.count),
    Ystd(re.Ystd),
    debug(re.debug)
    {   
    }

    // // operator=
    ReactionEntry& operator= (const ReactionEntry& re)
    {
        Y = re.Y;
        T = re.T;
        p = re.p;
        dtChem = re.dtChem;
        rho0 = re.rho0;
        rho = re.rho;
        count = re.count;
        Ystd = re.Ystd;
        debug = re.debug;
        return *this;
    }

    //- data member
    scalarField Y;
    scalar      T;
    scalar      p;
    scalar      dtChem;
    scalar      rho0;
    scalar      rho;
    label       count;
    scalarField Ystd;
    bool        debug;

    void read(Istream& is)
    {
        is >> Y >> T >> p >> dtChem >> rho0 >> rho >> count >> Ystd >> debug;
    }

    void write(Ostream& os) const
    {
        os << Y << " ";
        os << T << " ";
        os << p << " ";
        os << dtChem << " ";
        os << rho0 << " ";
        os << rho << " ";
        os << count << " ";
        os << Ystd << " ";
        os << debug << endl;
    }

};

}
#endif

