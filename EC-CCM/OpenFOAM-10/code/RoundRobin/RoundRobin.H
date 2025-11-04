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

#ifndef ROUND_ROBIN_H
#define ROUND_ROBIN_H

#include "fvCFD.H"

namespace Foam
{

// declaration


class RoundRobin;


class RoundRobin
{

    public:
    // // constructor
    //- Construct from the number of cores, the LCS include all cores.
    RoundRobin
    (
        const label& nCores
    ):
    nCores_(nCores),
    myCoreNum_(Pstream::myProcNo()),
    myCoreNumLCS_(Pstream::myProcNo()),
    myPosition_(Pstream::myProcNo()),
    coreList_(labelList(nCores,0)),
    teamsOrder_(labelList(nCores,0))
    {
        initialization();
        initializeDummyCore();
    }

    // Constructor LCS every N cores
    RoundRobin
    (
        const label& nCores,
        const label& nCoresLCS
    ):
    nCores_(nCoresLCS),
    myCoreNum_(0),
    myCoreNumLCS_(0),
    myPosition_(0),
    coreList_(labelList(nCoresLCS,0)),
    teamsOrder_(labelList(nCoresLCS,0))
    {

        // check nCores is a multiple of nCoresLCS
        if (nCores % nCoresLCS != 0)
        {
            FatalErrorIn("RoundRobin::RoundRobin")
            << "Number of cores must be a multiple of the LCS cores"
            << abort(FatalError);
        }
        initialization();
        

        // LCS != GCS, treat coreList_
        myCoreNum_ = Pstream::myProcNo();
        myPosition_ = myCoreNum_ % nCoresLCS;
        myCoreNumLCS_ = myPosition_;
        
        label round = (myCoreNum_-myPosition_)/nCoresLCS;
        forAll (coreList_, i)
        {
            coreList_[i] = round*nCoresLCS + i;
        }
        initializeDummyCore();



        //Pout << "coreList_:" << coreList_ << endl;
        
        
        
    }


    void initializeDummyCore()
    {
        if (oddCores_)
        {
            coreList_[nEffCores_-1] = -1;
        }
    }

    //- initialize oddCores_/nEffCores_/teamsOrder_/opponents
    //- initialize coreList_ correctly if LCS is the GCS
    //- If LCS != GCS, coreList_ should be treated properly
    void initialization()
    {
        if (nCores_ < 2)
        {
            FatalErrorIn("RoundRobin::RoundRobin")
            << "Number of cores must be greater than 1"
            << abort(FatalError);
        }


        oddCores_ = (nCores_ % 2 == 0) ? false : true;
        nEffCores_ = (nCores_ % 2 == 0) ? nCores_ : nCores_ + 1;
        if (oddCores_)
        {
            coreList_.resize(nEffCores_);
            teamsOrder_.resize(nEffCores_);
        }
        forAll(coreList_,i)
        {
            coreList_[i] = i;
            teamsOrder_[i] = i;
        }

        teamsOrder_ = coreList_;
        opponents_ = labelList(totalRounds(),0);
        opponentsLCS_ = labelList(totalRounds(),0);
    }
    
    
    void shift()
    {
        // except from the first element, all the other elements are shifted to the right
        // the last one goes to the second
        labelList newTeamsOrder(nEffCores(),0);
        newTeamsOrder[0] = teamsOrder_[0];
        newTeamsOrder[1] = teamsOrder_[nEffCores_-1];
        forAll(teamsOrder_,i)
        {
            if (i > 1)
            {
                newTeamsOrder[i] = teamsOrder_[i-1];
            }
        }
        teamsOrder_ = newTeamsOrder;

        if (myPosition_!=0)
        {
            myPosition_ ++;
            if (myPosition_ == nEffCores())
            {
                myPosition_ = 1;
            }
        }
    }

    label competeWith(label coreNo)
    {
        return nEffCores_ - coreNo - 1;
    }

    label nEffCores() const
    {
        return nEffCores_;
    }

    labelList teamsOrder() const
    {
        return teamsOrder_;
    }

    label myCoreNum() const
    {
        return myCoreNum_;
    }

    label totalRounds() const
    {
        return nEffCores()-1;
    }

    void recordOpponents(bool debug = false)
    {
        for (label i=0; i!=totalRounds()-1;i++)
        {
            opponents_[i] = coreList_[teamsOrder_[competeWith(myPosition_)]];
            opponentsLCS_[i] = teamsOrder_[competeWith(myPosition_)];
            shift();

            //Info << "teamsOrder_:" << teamsOrder_ << endl; 

        }
        //Info << "teamsOrder_:" << teamsOrder_ << endl;
        opponents_[totalRounds()-1] = coreList_[teamsOrder_[competeWith(myPosition_)]];
        opponentsLCS_[totalRounds()-1] = teamsOrder_[competeWith(myPosition_)];
        if (debug)
        {
            Pout << "Opponents for core " << myCoreNum_ << " are: " << opponents_ << endl;
            Pout << "Opponents in LCS for core " << myCoreNum_ << " are: " << opponentsLCS_ << endl;
            Pout << "myCoreNum_:" << myCoreNum_ << endl;
            Pout << "myCoreNumLCS_:" << myCoreNumLCS_ << endl;
        }
    }

    // //  default constructor
    RoundRobin():
    RoundRobin(2)
    {}




    // // copy constructor


    // // data member

    //- The number of cores of the local communication system (LCS)
    label nCores_;

    //- The core number of the current processor in GCS
    label myCoreNum_;

    //- The core number in LCS
    label myCoreNumLCS_;

    //- Current position in LCS
    label myPosition_;

    //- All core numbers in the global communication system (GCS)
    labelList coreList_;

    //- The order of teams in the LCS
    labelList teamsOrder_;

    //- Whether the number of cores is odd
    bool oddCores_;

    //- Effective cores, if the number of cores is odd, it is nCores_+1
    //- In the odd case, the last core is dummy
    label nEffCores_;

    //- Opponents for each round for the current core
    labelList opponents_;

    //- Opponents in LCS
    labelList opponentsLCS_;



};

}
#endif

