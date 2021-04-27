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
    solid4GeneralContactFvPatchVectorField

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
    solid4GeneralContactFvPatchVectorField.C
    solid4GeneralContactFvPatchVectorFieldCalcContact.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Rewrite by Hrvoje Jasak.

\*---------------------------------------------------------------------------*/

#ifndef solid4GeneralContactFvPatchVectorField_H
#define solid4GeneralContactFvPatchVectorField_H

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
               Class solid4GeneralContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solid4GeneralContactFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    // Private data
		
	//*************** based on solidGeneral*****************
		//- Is this the global master patch: the global master patches keeps a
        // copy of mesh.allPoints() and updates them to the deformed position
        // each time before the contact is corrected
        mutable bool* firstPatchPtr_;
		
		//- Patch index of global master patch
        mutable label* firstPatchIndexInListPtr_;
		
		//- Is this the current master patch
        mutable bool* currentMasterPtr_;
		
		//- Patch index of current master patch
        mutable label* currentMasterIndexPtr_;
		
		//- If the current patch is the local slave of the contact pair
        // The local slave is the patch with the lower patch index
        mutable boolList* localSlavePtr_;
		
		//- List of slaveG patch names
        mutable List<word>* slaveGPatchNamesPtr_;

        //YYY- List of slaveG patch indices 
        mutable List<label>* slaveGPatchIndicesPtr_;
		
		//- List of slaveG zone names
        mutable List<word>* slaveGZoneNamesPtr_;

        //- List of slaveG zone indices
        mutable List<label>* slaveGZoneIndicesPtr_;
	//*************** END based on solidGeneral**************
	
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
    //    mutable wordList slavePatchNames_;
		
		//************ based on solidGeneralContact****************
		mutable List<word>* slavePatchNames_;
		//************ END based on solidGeneralContact****************

        //- Slave patch indices
        mutable labelList* slavePatchIndicesPtr_;

        //- Is the contact boundary active
        Switch rigidMaster_;

        //- Normal contact model pointers
        mutable PtrList<generalNormalContactModel> normalModels_;

        //- Friction contact model pointer
        mutable PtrList<generalFrictionContactModel> frictionModels_;

        //- Store current penalty scale factors (if using penalty contact model)
        mutable scalarField normalPenaltyFactors_;

        //- Master face zone kept in deformed configuration
        mutable globalPolyPatch* zonePtr_;

        //- Slave face zones kept in deformed configuration
        //  Stored by the master patch
        mutable PtrList<globalPolyPatch> slaveZones_;

        //- Zone-to-zone interpolations
        mutable PtrList<newGgiStandAlonePatchInterpolation> zoneToZones_;

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
		
		//*************** based on solidGeneral*****************
		//- Offset used if bounding box of patch has a zero dimension
        mutable scalar bbOffset_;
		//*************** END based on solidGeneral*****************


    // Private Member Functions

        //- Is a moving mesh (e.g. updated Lagrangian) approach be used
        bool movingMesh() const;

        //- Make the slave patch names
        //void makeSlavePatchNames(const dictionary& dict) const;
		
		// ************ based on solidGeneral ****New declaration ***
		//- Make the slave patch names
        void makeSlavePatchNames() const;

        //- Calculate the slave patch indices
        void calcSlavePatchIndices() const;

        //- Make normal contact models
        void makeNormalModels(const dictionary& dict) const;

        //- Make friction contact models
        void makeFrictionModels(const dictionary& dict) const;

        //- Clear out demand driven data
        //virtual void clearOut();

        //- Create the zone
        void calcZone() const;

        //- Create the slave zones
        void calcSlaveZones() const;

        //- Create the zone-to-zone interpolators
        void calcZoneToZones() const;

        //- Move the master and slave zones to the deformed configuration
        void moveZonesToDeformedConfiguration();
		
        // Set the contactPerSlave field
        void calcContactPerSlave() const;
		
		//*************** based on solidGeneral*****************
		//- Check global master patch
        void calcFirstPatchInList() const;

        //- Check global master patch index
        void calcFirstPatchIndexInList() const;
		
		//- Check current master patch
        void calcCurrentMaster() const;
		
		//- Check current master patch index
        void calcCurrentMasterIndex() const;
		
		//- Calculate local slaves
        void calcLocalSlave() const;
		
		//- Return local slaves
			const boolList& localSlave() const;
			
		//- Calculate slaveG patch names
        void calcslaveGPatchNames() const;

        //- Calculate slaveG zone names
        void calcslaveGZoneNames() const;
		//*************** END based on solidGeneral*****************

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
		
		//*************** based on solidGeneral*****************
		//- Calculate the bounding box offset
        void calcBbOffset() const;

        //- Return the bounding box offset
        scalar bbOffset() const;
		//*************** END based on solidGeneral*****************


public:

    //- Runtime type information
    TypeName("solid4GeneralContact");


    // Constructors

        //- Construct from patch and internal field
        solid4GeneralContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solid4GeneralContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solid4GeneralContactFvPatchVectorField onto
        //  a new patch
        solid4GeneralContactFvPatchVectorField
        (
            const solid4GeneralContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solid4GeneralContactFvPatchVectorField
        (
            const solid4GeneralContactFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solid4GeneralContactFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solid4GeneralContactFvPatchVectorField
        (
            const solid4GeneralContactFvPatchVectorField&,
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
                new solid4GeneralContactFvPatchVectorField(*this, iF)
            );
        }

    //- Destructor
    virtual ~solid4GeneralContactFvPatchVectorField();


    // Member functions

        // Access
			
			//*************** based on solidGeneral*****************
			//- Return true if this patch is the global master
            bool firstPatchInList() const;

            //- Return index of global master patch
            label firstPatchIndexInList() const;
			
			//- Return true if this patch is the current master
            bool currentMaster() const;

            //- Return index of current master patch
            label currentMasterIndex() const;			
			
			//Step1 - Return the list of slaveG patch names
            const List<word>& slaveGPatchNames() const;

            //Step3 - Return the list of slaveG patch indices
            const List<label>& slaveGPatchIndices() const;
        
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
			//*************** END based on solidGeneral*****************
			
            //- Return if this patch is the slave
            Switch master() const
            {
                return master_;
            }

            //- Return the slave patch names
            //const wordList& slavePatchNames() const;
			
			//*************** based on solidGeneral*****************
			//- Return the list of slave patch names
            const List<word>& slavePatchNames() const;
			//*************** END based on solidGeneral*****************

            //- Return the slave patch indices
            const labelList& slavePatchIndices() const;

            //- Return the slave patch
            //  This function is not allowed if slavePatchIndices.size() > 1
            const solid4GeneralContactFvPatchVectorField& slavePatchField() const;

            //- Return reference to the normal contact models
            PtrList<generalNormalContactModel>& normalModels();

            //- Return const reference to the normal contact models
            const PtrList<generalNormalContactModel>& normalModels() const;

            //- Return reference to the friction contact models
            PtrList<generalFrictionContactModel>& frictionModels();

            //- Return const reference to the friction contact models
            const PtrList<generalFrictionContactModel>& frictionModels() const;

            //- Return normal contact model for the slave patch calling this
            //  function
            generalNormalContactModel& normalModelForThisSlave();

            //- Return friction contact model for the slave patch calling this
            //  function
            generalFrictionContactModel& frictionModelForThisSlave();

            //- Return const reference to deformed zone
            const globalPolyPatch& zone() const;

            //- Return reference to deformed zone
            globalPolyPatch& zone();

            //- Return const reference to deformed slave zones
            const PtrList<globalPolyPatch>& slaveZones() const;

            //- Return reference to deformed slave zones
            PtrList<globalPolyPatch>& slaveZones();

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
			
			//************** based on solidGeneral************
			//- Given the patch index of a slave, return the position in the
            // slave list of the current patch
            label findSlaveID(const label patchID) const;
			//************** END based on solidGeneral**********

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

            //- Clear out demand driven data
            virtual void clearOut();

            //- Return the rate of dissipated energy due to friction for the
            //  master patch
            virtual tmp<scalarField> frictionHeatRate() const;

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

//#ifdef NoRepository
//#   include "solid4GeneralContactFvPatchVectorFieldTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
