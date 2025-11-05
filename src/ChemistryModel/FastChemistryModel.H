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

Class
    Foam::FastChemistryModel

Description

SourceFiles
    FastChemistryModelI.H
    FastChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef FastChemistryModel_H
#define FastChemistryModel_H

#include <chrono>

#include "basicFastChemistryModel.H"
#include "dataBlock.H"
#include "simpleDataBlock.H"
#include "OptReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FastChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class FastChemistryModel
:
    public basicFastChemistryModel
{
    // Private data

        // //- Reference to the field of specie mass fractions
        // const PtrList<volScalarField>& Yvf_;

        //- Number of species
        label nSpecie_;

        //- Type of the Jacobian to be calculated
        const jacobianType jacobianType_;

        //- Reactions
        mutable OptReaction reaction;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField> RR_;

public:
        //- Data buffer
        mutable double*     buffer;
        mutable double*     YTpWork[12];
        mutable double*     YTpYTpWork[3];
        
        //- Size of the ODE system
        label n_;
        mutable unsigned int alignN=0;

private:

        //- Minimum reaction temperature
        const scalar Treact;

            //****************************MPI Communication*****************************//

            //- CPU time larger tave*DLBthreshold is busy process, and will send data to other process.
            scalar DLBthreshold = 1.0;

            //- Updte load balancing per Iter
            label MaxIter = 1;
            label CurrentIter = 0;

            //- CPU time for data transfer
            List<List<int64_t>> cpuLoadTransferTable;

            //- CPU time to solve chemistry for each cell
            std::vector<std::pair<int64_t,label>> CPUtimeField;

            //- CPU time so solve chemistry for each process
            List<int64_t> chemistryIntegrationTime;

            //- Size of send buffer
            List<std::streamsize> sendBufferSize_;

            //- Size of receive buffer
            List<std::streamsize> recvBufferSize_;

            //- Send buffer, containing the binary data block
            List<DynamicList<char>> sendBuffer_;

            //- Receive buffer, containing the binary data block
            List<DynamicList<char>> recvBuffer_;

            //- The receive position of receive buffer
            List<label> recvBufPos_;

            //- firt time to be executed
            bool firstTime = true;

            //- Skip this cell or solve?
            List<bool> skip;

            //- Determine whether this process is busy
            List<bool> IamBusyProcess;

            //- Whether to switch on the load balancing
            bool Balance;

            void exchange
            (
                const UList<DynamicList<char>>& sendBufs,
                const List<std::streamsize>& recvSizes,
                List<DynamicList<char>>& recvBufs,
                const int tag,
                const label comm,
                const bool block
            );

            void exchangeSizes
            (
                const UList<DynamicList<char>>& sendBufs,
                labelList& recvSizes,
                const label comm
            );
            //****************************MPI Communication*****************************//

    // Private Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField>& RR();

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT){return 0;}

public:

    //- Runtime type information
    TypeName("FastChemistryModel");


    // Constructors

        //- Construct from mesh
        FastChemistryModel(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        FastChemistryModel(const FastChemistryModel&) = delete;


    //- Destructor
    virtual ~FastChemistryModel();


    // Member Functions

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Calculates the reaction rates
        virtual void calculate();

        // Return the reaction object
        inline OptReaction&  getReaction() const;

        // Chemistry model functions (overriding abstract functions in
        // basicFastChemistryModel.H)

            //- The number of species
            inline virtual label nSpecie() const;

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m^3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );

            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT) ;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT) ;

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;

        // ODE functions (overriding abstract functions in ODE.H)
            virtual void derivatives
            (
                const scalar t,
                const label li,
                const double p,
                double* __restrict__ Phi,
                double* __restrict__ dPhidt,
                double* __restrict__ Cp,
                double* __restrict__ Ha
            ) const ;

            virtual void jacobian
            (
                const scalar t,
                const label li,
                const double p,
                double* __restrict__ Phi,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac
            ) const ;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& Y,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;

            virtual void solve
            (
                const label li,
                double T,
                double& __restrict__ deltaT,
                double& __restrict__ subDeltaT
            ) const = 0;  

        // Rate evaluation on the basis of cells
            scalarField getRRGivenYTP
            (
                const scalarField& Y,
                const scalar T,
                const scalar p,
                const scalar deltaT,
                scalar& deltaTChem,
                const scalar& rho,
                const scalar& rho0
            ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FastChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FastChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
