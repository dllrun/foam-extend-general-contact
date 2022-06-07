/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    solidGeneralContactFvPatchVectorField

Description
    Contact boundary conditions for finite volume solidMechanics solvers.

    The treatment of the normal contact and tangential contact can be specified
    individually at run-time using the generalNormalContactModel and
    generalFrictionContactModel.

    The boundary condition keeps a copy of the master and slave patch and
    keeps them in the deformed configuration; this allows the calculation of
    contact distances and interpolation between the globalPatchZones.

    The distance calculations and interpolations are performed by the GGI class.

    More details in:

    P. Cardiff, A. Karać, A. Ivanković: Development of a Finite Volume contact
    solver based on the penalty method. Computational Materials Science,
    03/2014, 64:283–284. DOI: 10.1016/j.commatsci.2012.03.011.

    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

    P. Cardiff, Development of the Finite Volume Method for Hip Joint Stress
    Analysis, PhD thesis, University College Dublin, 2012.


SourceFiles
    solidGeneralContactFvPatchVectorField.C

Author
	Ranjan Dhakal, IPPE, Graz University of Technology 
Co-Author	 
	Stefan Radl, IPPE, Graz University of Technology  
    Philip Cardiff, University College Dublin 

\*---------------------------------------------------------------------------*/

#ifndef solidGeneralContactFvPatchVectorField_H
#define solidGeneralContactFvPatchVectorField_H

#include "foamTime.H"
#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "generalNormalContactModel.H"
#include "generalFrictionContactModel.H"
#include "globalPolyPatch.H"
#include "newGgiInterpolation.H"
#include "Switch.H"
#include "boundBox.H"
#include "standAlonePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class solidGeneralContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidGeneralContactFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    // Private data
		

		//- Is this the global master patch: the global master patches keeps a
        // copy of mesh.allPoints() and updates them to the deformed position
        // each time before the contact is corrected
        mutable bool* firstPatchInListPtr_;
		
		//- Patch index of global master patch
        mutable label* firstPatchInListIndexPtr_;
		
		//- Is this the current master patch
        mutable bool* currentMasterPtr_;
		
		//- Patch index of current master patch
        mutable label* currentMasterIndexPtr_;
		
		//- Is this the current slave patch
        mutable bool* currentSlavePtr_;
		
		//- If the current patch is the local slave of the contact pair
        // The local slave is the patch with the lower patch index
        mutable boolList* localMasterPtr_;
		
		//- List of slaveG patch names
        mutable List<word>* slaveGPatchNamesPtr_;

        //YYY- List of slaveG patch indices 
        mutable List<label>* slaveGPatchIndicesPtr_;
		
		//- List of slaveG zone names
        mutable List<word>* slaveGZoneNamesPtr_;

        //- List of slaveG zone indices
        mutable List<label>* slaveGZoneIndicesPtr_;

	
        //- Store a copy of the patch dictionary
        //  This allows us to use lazy evaluation when creating the contact
        //  models
        const dictionary dict_;

        //- Is this patch the master
        const Switch master_;

        //- Write the deformed contact zones to VTK surface files
        const Switch writeZoneVTK_;

        //- Write the contact point distances fields for viewing
        const Switch writePointDistanceFields_;

        //- Slave patch names
    //    mutable wordList targetPatchNamesPtr_;
		
		mutable List<word>* targetPatchNamesPtr_;


        //- Slave patch indices
        mutable labelList* targetPatchIndicesPtr_;

        //- Is the contact boundary active
        Switch rigidMaster_;

        //- Normal contact model pointers
        mutable PtrList<generalNormalContactModel> normalModels_;

        //- Friction contact model pointer
        mutable PtrList<generalFrictionContactModel> frictionModels_;

        //- Store current penalty scale factors (if using penalty contact model)
        mutable scalarField normalPenaltyFactors_;

        //- Master face zone kept in deformed configuration
        mutable globalPolyPatch* globalPolyPatchFaceZonePtr_;

        //- Slave face zones kept in deformed configuration
        //  Stored by the master patch
        mutable PtrList<globalPolyPatch> targetGlobalPolyPatchFaceZones_;

        //- Zone-to-zone interpolations
        mutable PtrList<newGgiStandAlonePatchInterpolation> zoneToZones_;
		

		//- Patch to patch algorithm
        const intersection::algorithm alg_;
		
		//- Patch to patch direction
        const intersection::direction dir_;


        //- GGI quick reject algorithm
        newGgiInterpolation::quickReject quickReject_;

        //- Top corner of the region where te contact is expected to be
        // This aims to save time in the GGI distance calculation
        const vector regionOfInterestTopCorner_;

        //- Bottom corner of the region where te contact is expected to be
        // This aims to save time in the GGI distance calculation
        const vector regionOfInterestBottomCorner_;

        //- Region where the contact is expected to be
        //  This aims to save time in the GGI distance calculations
        const boundBox regionOfInterest_;

        //- Contact field for all the contacts, stored on the current patch
        //  1 means in contact
        //  0 means not in contact
        scalarField contact_;

        //- Contact field corresponding to each slave patch, stored on the
        //  current patch
        //  1 means in contact
        //  0 means not in contact
        mutable PtrList<scalarField> contactPerSlave_;

        //- Flag to enable approach for scaling traction near the downstream
        //  patch to fix an issue with unphysical behaviour
        Switch scaleFaceTractionsNearDownstreamPatch_;

        //- Field to scale traction of slave patch
        mutable autoPtr<scalarField> scaleTractionFieldPtr_;

        //- Current time index
        label curTimeIndex_;
		

		//- The current patch tractions for each contatc pair
        mutable List<vectorField>* curSlaveTractionPtr_;
		

		//- The current patch tractions for each contatc pair
        mutable List<vectorField>* curPatchTractionPtr_;
		
		// Rate of friction energy dissipation
        mutable scalarField* QcPtr_;
		
		// Rate of friction energy dissipation for each contact pair
        mutable List<scalarField>* QcsPtr_;
		
		//- Offset used if bounding box of patch has a zero dimension
        mutable scalar bbOffset_;
		
		mutable scalarField slaveBbOffsetList_;

		
		// //- Zone index i.e. index of global face zone corresponding to the
        // current patch 
		mutable label zoneIndex_;
		
		//- List of slave zone names
        mutable List<word>* targetGlobalPolyPatchFaceZoneNamesPtr_;

        //- List of slave zone indices
        mutable List<label>* targetGlobalPolyPatchFaceZoneIndicesPtr_;		
			

    // Private Member Functions

        //- Is a moving mesh (e.g. updated Lagrangian) approach be used
        bool movingMesh() const;

        //- Make the slave patch names
        //void maketargetPatchNames(const dictionary& dict) const;
		
		// ************ based on solidGeneral ****New declaration ***
		//- Make the slave patch names
        void maketargetPatchNames() const;

        //- Calculate the slave patch indices
        void calcTargetPatchIndices() const;
		
		
		void calcNormalContactModels(const dictionary& dict) const;
		
		

        //- Make normal contact models
        void calcNormalModels(const dictionary& dict) const;

        //- Make friction contact models
        void calcFrictionModels(const dictionary& dict) const;

        //- Clear out demand driven data
        //virtual void clearOut();


		//- Calculate slave patch names
        void calcTargetPatchNames() const;

        //- Calculate slave zone names
        void calcTargetGlobalPolyPatchFaceZoneNames() const;
		
		//- Calculate zone index
        void calcGlobalPolyPatchFaceZoneIndex() const;


        //- Create the zone
        void calcGlobalPolyPatchFaceZone() const;

        //- Create the slave zones
        void calcTargetGlobalPolyPatchFaceZones() const;

        //- Create the zone-to-zone interpolators
        void calcGlobalPolyPatchFaceZoneToZones() const;
		
		//- Create the zone-to-zone interpolators
        void makeZoneToZone() const;

        //- Move the master and slave zones to the deformed configuration
        void moveZonesToDeformedConfiguration();
		
        // Set the contactPerSlave field
        void calcContactPerSlave() const;
		

		//- Check global master patch
        void calcFirstPatchInList() const;

        //- Check global master patch index
        void calcFirstPatchInListIndex()const;
		
		//- Check current master patch
        void calcCurrentMaster() const;
		
		//- Check current master patch index
        void calcCurrentMasterIndex() const;
		
		//- Check current slave patch
        void calcCurrentSlave() const;
		
		//- Calculate local slaves
        void calcLocalSlave() const;
		
		//- Return local slaves
			const boolList& localMaster() const;
			
		//- Calculate slaveG patch names
        void calcslaveGPatchNames() const;

        //- Calculate slaveG zone names
        void calcslaveGZoneNames() const;
		
		//- Make the current patch traction list
        void makeCurPatchTractions() const;
		
		//- Calculate rate of dissipated friction energy
        void calcQc() const;
		
		//- Calculate rate of dissipated friction energy for each contact
        //  pair
        void calcQcs() const;

		

        // //- Interpolate zone face values to points
        // template<class Type>
        // tmp<Field<Type> > zoneFaceToPointInterpolate
        // (
        //     const label zoneIndex,
        //     const Field<Type>& zoneField
        // ) const;

        // //- Interpolate zone point values to faces
        // template<class Type>
        // tmp<Field<Type> > zonePointToFaceInterpolate
        // (
        //     const label zoneIndex,
        //     const Field<Type>& zonePointField
        // ) const;

        //- Return field to scale the traction field on the slave patch
        scalarField scaleTractionField() const;

        //- Make the scaleTractionField
        void makeScaleTractionField() const;
		

		//- Calculate the bounding box offset
        void calcBbOffset() const;
		
		void calcSlaveBbOffset() const;

        //- Return the bounding box offset
        scalar bbOffset() const;
		
		scalarField slaveBbOffset() const;
		
		
		//- Get zone field from patch field
        template<class Type>
        tmp<Field<Type> > zoneField
        (
            const label zoneIndex,
            const label patchIndex,
            const Field<Type>& patchField
        ) const;
		
		
		//- Get patch field from zone field
        template<class Type>
        tmp<Field<Type> > patchField
        (
            const label patchIndex,
            const label zoneIndex,
            const Field<Type>& zoneField
        ) const;
		
		//- Interpolate zone face values to points
        template<class Type>
        tmp<Field<Type> > zoneFaceToPointInterpolate
        (
            const label zoneIndex,
            const Field<Type>& zoneField,
            const label shadowI
        ) const;

        //- Interpolate zone point values to faces
        template<class Type>
        tmp<Field<Type> > zonePointToFaceInterpolate
        (
            const label zoneIndex,
            const Field<Type>& zonePointField,
            const label shadowI
        ) const;
		


public:

    //- Runtime type information
    TypeName("solidGeneralContact");


    // Constructors

        //- Construct from patch and internal field
        solidGeneralContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidGeneralContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidGeneralContactFvPatchVectorField onto
        //  a new patch
        solidGeneralContactFvPatchVectorField
        (
            const solidGeneralContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidGeneralContactFvPatchVectorField
        (
            const solidGeneralContactFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidGeneralContactFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidGeneralContactFvPatchVectorField
        (
            const solidGeneralContactFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidGeneralContactFvPatchVectorField(*this, iF)
            );
        }

    //- Destructor
    virtual ~solidGeneralContactFvPatchVectorField();


    // Member functions

        // Access
		
		//- Return local slaves
			const boolList& localTSlave() const
			{
                return localMaster();
            }

			
			//- Return true if this patch is the global master
            bool firstPatchInList() const;

            //- Return index of global master patch
            label firstPatchInListIndex() const;
			
			//- Return true if this patch is the current master
            bool currentMaster() const;

            //- Return index of current master patch
            label currentMasterIndex() const;

			//- Return true if this patch is the current slave
            bool currentSlave() const;
			
			//Step1 - Return the list of slaveG patch names
            const List<word>& slaveGPatchNames() const;

            //Step3 - Return the list of slaveG patch indices
            const List<label>& slaveGPatchIndices() const;


			//- Return index of the zone
            label globalPolyPatchFaceZoneIndex() const;
			
			//- Return the list of slave zone names
            const List<word>& targetGlobalPolyPatchFaceZoneNames() const;

            //- Return the list of slave zone indices
            const List<label>& targetGlobalPolyPatchFaceZoneIndices() const;

        
            //- Return the list of slaveG zone names
            const List<word>& slaveGZoneNames() const;

            //- Return the list of slaveG zone indices
            const List<label>& slaveGZoneIndices() const;
			
			//Step2's - Return const reference to deformed master zone
            const standAlonePatch& slaveGZone(const label slaveGI) const;

            //- Return reference to deformed master zone
            standAlonePatch& slaveGZone(const label slaveGI);
			
			//- Given the patch index of a slaveG, return the position in the
            // slaveG list of the current patch
            label findslaveGID(const label patchID) const;
			
			//- Return the current patch traction for slaveI
            const vectorField& curPatchTractions(const label slaveI) const;
			
			//- Return the current patch traction for slaveI
            vectorField& curPatchTractions(const label slaveI);
			
			
			const vectorField& curSlaveTractions(const label slaveI) const;
			
			//- Return the current patch traction for slaveI
            vectorField& curSlaveTractions(const label slaveI);
			
			
            //- Return if this patch is the slave
            Switch master() const
            {
                return master_;
            }

            //- Return the slave patch names
            //const wordList& targetPatchNames() const;
			
			//- Return the list of slave patch names
            const List<word>& targetPatchNames() const;

            //- Return the slave patch indices
            const labelList& targetPatchIndices() const;

            //- Return the slave patch
            //  This function is not allowed if shadowPatchIndices.size() > 1
            const solidGeneralContactFvPatchVectorField& shadowPatchField() const;


            //- Return reference to the normal contact model
            //PtrList<generalNormalContactModel>& 
			generalNormalContactModel& normalContactModel(const label shadowI);

            //- Return const reference to the normal contact model
            //PtrList<generalNormalContactModel>& 
			const generalNormalContactModel& normalContactModel(const label shadowI) const;

			

			//- Return reference to the normal contact models
            generalNormalContactModel& normalModel(const label shadowI);

            //- Return const reference to the normal contact models
            const generalNormalContactModel& normalModel(const label shadowI) const;

            //- Return reference to the friction contact models
            generalFrictionContactModel& frictionModel(const label shadowI);

            //- Return const reference to the friction contact models
            const generalFrictionContactModel& frictionModel(const label shadowI) const;

            //- Return normal contact model for the slave patch calling this
            //  function
            generalNormalContactModel& normalModelForThisSlave(const label shadowI);

            //- Return friction contact model for the slave patch calling this
            //  function
            generalFrictionContactModel& frictionModelForThisSlave(const label shadowI);
			
			
            //- Return const reference to deformed zone
            const globalPolyPatch& globalPolyPatchFaceZone() const;

            //- Return reference to deformed zone
            globalPolyPatch& globalPolyPatchFaceZone();
			
			
			//- Return const reference to deformed master zone
            const globalPolyPatch& targetGlobalPolyPatchFaceZone(const label shadowI) const;

            //- Return reference to deformed master zone 
            globalPolyPatch& targetGlobalPolyPatchFaceZone(const label shadowI);
			
			//- Return const reference to patch-to-patch interpolation
            const newGgiStandAlonePatchInterpolation& zoneToZone
            (
                const label shadowI
            ) const;

            //- Return reference to patch-to-patch interpolation
            newGgiStandAlonePatchInterpolation& zoneToZone
            (
                const label shadowI
            );	
			
			

            //- Return const reference to deformed slave zones
            const PtrList<globalPolyPatch>& targetGlobalPolyPatchFaceZones() const;

            //- Return reference to deformed slave zones
            PtrList<globalPolyPatch>& targetGlobalPolyPatchFaceZones();

            //- Return const reference to patch-to-patch interpolation
            const PtrList<newGgiStandAlonePatchInterpolation>&
                zoneToZones() const;

            //- Return reference to patch-to-patch interpolation
            PtrList<newGgiStandAlonePatchInterpolation>& zoneToZones();

            //- Return const reference to zone-to-zone interpolation for the
            //  current slave patch; an error is thrown if the master calls this
            //  function
            const newGgiStandAlonePatchInterpolation&
            zoneToZoneForThisSlave() const;

            //- Return the global patch for the current slave patch
            //  An error is thrown if the master calls this function
            const globalPolyPatch& zoneForThisSlave() const;

            //- Return reference to contact field
            const scalarField& contact() const
            {
                return contact_;
            }

            //- Return reference to contactPerSlave fields
            PtrList<scalarField>& contactPerSlave();

            //- Return const reference to contactPerSlave fields
            const PtrList<scalarField>& contactPerSlave() const;
			
			
			//- Given the patch index of a slave, return the position in the
            // slave list of the current patch
            label findShadowID(const label patchID) const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // Member functions
		
			//- Dissipated energy due to friction
            virtual const scalarField& Qc() const;

            //- Dissipated energy due to friction from contact with shadowI
            virtual const scalarField& Qc(const label shadowI) const;
			

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

//#ifdef NoRepository
//#   include "solidGeneralContactFvPatchVectorFieldTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
