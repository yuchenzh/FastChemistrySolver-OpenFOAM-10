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

#ifndef RR_ENTRY_H
#define RR_ENTRY_H

#include "fvCFD.H"

namespace Foam
{

// declaration


class RREntry;
Ostream& operator << (Ostream& os, const RREntry& re);

Istream& operator >> (Istream& is, RREntry& re);

bool operator == (const RREntry& re1, const RREntry& re2);

bool operator != (const RREntry& re1, const RREntry& re2);





class RREntry
{

    public:
    //- constructor
    RREntry
    (
        const scalarField& RRvalue,
        const scalar& deltaTChemValue,
        const bool debugValue = false,
        const scalarField& Yvalue = scalarField(0)
    ):
    RR(RRvalue),
    deltaTChem(deltaTChemValue),
    debug(debugValue),
    Y(Yvalue)
    {}

    // // construct from Foam::Istream&
    RREntry(Istream& is)
    {
        read(is);
    }

    // //  default constructor
    RREntry():
    RR(scalarField()),
    deltaTChem(0),
    debug(false),
    Y(scalarField(0))
    {}

    // // copy constructor
    RREntry(const RREntry& re):
    RR(re.RR),
    deltaTChem(re.deltaTChem),
    debug(re.debug),
    Y(re.Y)
    {}

    // // operator=
    RREntry& operator= (const RREntry& re)
    {
        RR = re.RR;
        deltaTChem = re.deltaTChem;
        debug = re.debug;
        Y = re.Y;
        return *this;
    }

    //- data member
    scalarField RR;
    scalar deltaTChem;
    bool   debug;
    scalarField Y;

    void read(Istream& is)
    {
        is >> RR >> deltaTChem >> debug >> Y;
    }

    void write(Ostream& os) const
    {
        os << RR << " ";
        os << deltaTChem << " ";
        os << debug << " ";
        os << Y << " " << endl;
    }
};

}
#endif

