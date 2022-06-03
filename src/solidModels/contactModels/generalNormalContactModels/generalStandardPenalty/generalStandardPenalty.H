/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    generalStandardPenalty

Description
    Standard penalty method for normal contact law.

    Using a fixed penalty factor, force is applied to the contact faces based on
    the gap.

    The penalty method is effectively like attaching springs to any face that
    comes into contact.

SourceFiles
    generalStandardPenalty.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef generalStandardPenalty_H
#define generalStandardPenalty_H

#include "generalNormalContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class generalStandardPenalty Declaration
\*---------------------------------------------------------------------------*/

class generalStandardPenalty
:
    public generalNormalContactModel
{
    // Private data

        //- Contact law dictionary
        mutable dictionary normalContactModelDict_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Slave pressure
        volVectorField slavePressureVolField_;

        //- Area in contact
        volScalarField areaInContactVolField_;

        //- Penalty factor
        mutable scalar penaltyFactor_;

        //- Penalty scale factor
        mutable scalar penaltyScale_;

        //- Under-relaxation factor
        const scalar relaxFac_;

        //- Current average penetration in the contact region
        //  This is zero if there is no contact
        scalar averagePenetration_;

        //- Current minimum (greatest) penetration in the contact region
        //  Remember that penetration is negative for points in contact
        scalar minPenetration_;

        //- Reference distance spring force calculation
        const scalar epsilon0_;

        //- Contact iteration number
        label contactIterNum_;


    // Private Member Functions

        //- Calculate penalty factor
        void calcPenaltyFactor() const;

        //- Disallow default bitwise assignment
        void operator=(const generalStandardPenalty&);


public:

    //- Runtime type information
    TypeName("generalStandardPenalty");


    // Constructors

        //- Construct from dictionary
        generalStandardPenalty
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const standAlonePatch& masterFaceZonePatch,
            const standAlonePatch& slaveFaceZonePatch
        );

        //- Construct as copy
        generalStandardPenalty(const generalStandardPenalty&);

        //- Construct and return a clone
        virtual autoPtr<generalNormalContactModel> clone() const
        {
            return autoPtr<generalNormalContactModel>(new generalStandardPenalty(*this));
        }


    // Destructor

        virtual ~generalStandardPenalty()
        {}


    // Member Functions


        //- Correct contact model
        virtual void correct
        (
            const vectorField& slavePatchFaceNormals,
            const scalarField& slavePointPenetration,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        );

        //- It is the start of a new time step e.g. update penalty factor
        virtual void newTimeStep() const
        {}

        //- Return slave pressure
        virtual const vectorField& slavePressure() const
        {
            return slavePressureVolField_.boundaryField()[slavePatchID()];
        }

        //- Return slave pressure
        virtual vectorField& slavePressure()
        {
			//Info<<"In -- generalStandardPenalty::slavePressure() line:"<<__LINE__<<endl;
            return slavePressureVolField_.boundaryField()[slavePatchID()];
        }

        //- Penalty factor
        virtual scalar penaltyFactor() const;

        //- Return area in contact of slave patch
        virtual const scalarField& areaInContact() const
        {
			Info<<"In generalStandardPenalty.H::areaInContact() line:"<<__LINE__<<endl;
			//Info<<"areaInContactVolField_ in generalStandardPenalty.H::areaInContact(): "<<areaInContactVolField_<<endl;
            Info<<"slavePatchID() in generalStandardPenalty.H::areaInContact(): "<<slavePatchID()<<endl;
			return areaInContactVolField_.boundaryField()[slavePatchID()];
			//***************** Testing master patchID ************************
            //return areaInContactVolField_.boundaryField()[masterPatchID()];
			//***************** End testing master patchID ************************
        };

        //- Return area in contact of slave patch
        virtual scalarField& areaInContact()
        {
			return areaInContactVolField_.boundaryField()[slavePatchID()];
			//***************** Testing master patchID ************************
            //return areaInContactVolField_.boundaryField()[masterPatchID()];
			//***************** End testing master patchID ************************
        };

        //- Update penalty scale factors (for penalty laws)
        virtual scalar updatePenaltyScale(const scalar previousPenaltyScale);

        //- Perform any reauired field mapping e.g. after a topo change
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
