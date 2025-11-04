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

#include "ParallelComm.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ParallelComm<ThermoType>::ParallelComm
(
    RoundRobin& RoundRobinCommunicator
)
:
    RoundRobinCommunicator_(RoundRobinCommunicator)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
template<typename DataType>
void Foam::ParallelComm<ThermoType>::sendReceivePair
(
    const label opponent,
    const label opponentLCS,
    const DataType& dataToSend,
    DataType& dataToReceive,
    const label myCoreNumLCS
)
{
    if (Pstream::parRun())
    {
        if (opponentLCS > myCoreNumLCS)
        {
            // Send first, then receive with proper scoping
            {
                OPstream sendMaterials(Pstream::commsTypes::scheduled, opponent);
                sendMaterials << dataToSend;
            }
            {
                IPstream receiveMaterials(Pstream::commsTypes::scheduled, opponent);
                receiveMaterials >> dataToReceive;
            }
        }
        else
        {
            // Receive first, then send with proper scoping
            {
                IPstream receiveMaterials(Pstream::commsTypes::scheduled, opponent);
                receiveMaterials >> dataToReceive;
            }
            {
                OPstream sendMaterials(Pstream::commsTypes::scheduled, opponent);
                sendMaterials << dataToSend;
            }
        }
    }
    else
    {
        // Serial run - just copy data
        dataToReceive = dataToSend;
    }
}


template<class ThermoType>
template<typename DataType>
void Foam::ParallelComm<ThermoType>::performRoundRobinDistribute
(
    const List<DataType>& sendData,
    List<DataType>& receiveData
)
{
    label myCoreNumLCS = RoundRobinCommunicator_.myCoreNumLCS_;
    
    // Send to myself (no network communication needed)
    receiveData[myCoreNumLCS] = sendData[myCoreNumLCS];
    
    // Communicate with other cores if parallel
    if (Pstream::parRun())
    {
        for (label i = 0; i < RoundRobinCommunicator_.totalRounds(); i++)
        {
            label opponentLCS = RoundRobinCommunicator_.opponentsLCS_[i];
            label opponent = RoundRobinCommunicator_.opponents_[i];
            bool dummyRound = (opponent == -1);

            if (!dummyRound)
            {
                sendReceivePair
                (
                    opponent,
                    opponentLCS, 
                    sendData[opponentLCS],
                    receiveData[opponentLCS],
                    myCoreNumLCS
                );
            }
        }
    }
    else
    {
        // Serial run - copy all data
        receiveData = sendData;
    }
}


template<class ThermoType>
template<typename DataType>
void Foam::ParallelComm<ThermoType>::performRoundRobinBroadcast
(
    DataType& myData,
    List<DataType>& allData,
    const label myCoreNumLCS
)
{
    // Initialize my own data in the list
    allData[myCoreNumLCS] = myData;
    
    // Communicate with other cores if parallel
    if (Pstream::parRun())
    {
        for (label i = 0; i < RoundRobinCommunicator_.totalRounds(); i++)
        {
            label opponentLCS = RoundRobinCommunicator_.opponentsLCS_[i];
            label opponent = RoundRobinCommunicator_.opponents_[i];
            bool dummyRound = (opponent == -1);

            if (!dummyRound)
            {
                sendReceivePair
                (
                    opponent,
                    opponentLCS,
                    myData,
                    allData[opponentLCS],
                    myCoreNumLCS
                );
            }
        }
    }
    else
    {
        // Serial run - all data is the same as my data
        forAll(allData, i)
        {
            allData[i] = myData;
        }
    }
}


// ************************************************************************* //