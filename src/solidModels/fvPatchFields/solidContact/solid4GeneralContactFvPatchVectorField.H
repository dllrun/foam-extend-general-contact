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
        mutable bool* globalMasterPtr_;
		
		//- Patch index of global master patch
        mutable label* globalMasterIndexPtr_;
		
		//- If the current patch is the local slave of the contact pair
        // The local slave is the patch with the lower patch index
        mutable boolList* localSlavePtr_;
		
		//- List of shadowG patch names
        mutable List<word>* shadowGPatchNamesPtr_;

        //YYY- List of shadowG patch indices 
        mutable List<label>* shadowGPatchIndicesPtr_;
		
		//- List of shadowG zone names
        mutable List<word>* shadowGZoneNamesPtr_;

        //- List of shadowG zone indices
        mutable List<label>* shadowGZoneIndicesPtr_;
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

        //- Shadow patch names
    //    mutable wordList shadowPatchNames_;
		
		//************ based on solidGeneralContact****************
		mutable List<word>* shadowPatchNames_;
		//************ END based on solidGeneralContact****************

        //- Shadow patch indices
        mutable labelList* shadowPatchIndicesPtr_;

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

        //- Shadow face zones kept in deformed configuration
        //  Stored by the master patch
        mutable PtrList<globalPolyPatch> shadowZones_;

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

        //- Contact field corresponding to each shadow patch, stored on the
        //  current patch
        //  1 means in contact
        //  0 means not in contact
        mutable PtrList<scalarField> contactPerShadow_;

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

        //- Make the shadow patch names
        //void makeShadowPatchNames(const dictionary& dict) const;
		
		// ************ based on solidGeneral ****New declaration ***
		//- Make the shadow patch names
        void makeShadowPatchNames() const;

        //- Calculate the shadow patch indices
        void calcShadowPatchIndices() const;

        //- Make normal contact models
        void makeNormalModels(const dictionary& dict) const;

        //- Make friction contact models
        void makeFrictionModels(const dictionary& dict) const;

        //- Clear out demand driven data
        //virtual void clearOut();

        //- Create the zone
        void calcZone() const;

        //- Create the shadow zones
        void calcShadowZones() const;

        //- Create the zone-to-zone interpolators
        void calcZoneToZones() const;

        //- Move the master and slave zones to the deformed configuration
        void moveZonesToDeformedConfiguration();
		
        // Set the contactPerShadow field
        void calcContactPerShadow() const;
		
		//*************** based on solidGeneral*****************
		//- Check global master patch
        void calcGlobalMaster() const;

        //- Check global master patch index
        void calcGlobalMasterIndex() const;
		
		//- Calculate local slaves
        void calcLocalSlave() const;
		
		//- Return local slaves
			const boolList& localSlave() const;
			
		//- Calculate shadowG patch names
        void calcshadowGPatchNames() const;

        //- Calculate shadowG zone names
        void calcshadowGZoneNames() const;
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
            bool globalMaster() const;

            //- Return index of global master patch
            label globalMasterIndex() const;
			
			//Step1 - Return the list of shadowG patch names
            const List<word>& shadowGPatchNames() const;

            //Step3 - Return the list of shadowG patch indices
            const List<label>& shadowGPatchIndices() const;
        
            //- Return the list of shadowG zone names
            const List<word>& shadowGZoneNames() const;

            //- Return the list of shadowG zone indices
            const List<label>& shadowGZoneIndices() const;
			
			//Step2's - Return const reference to deformed master zone
            const standAlonePatch& shadowGZone(const label shadowGI) const;

            //- Return reference to deformed master zone
            standAlonePatch& shadowGZone(const label shadowGI);
			
			//- Given the patch index of a shadowG, return the position in the
            // shadowG list of the current patch
            label findshadowGID(const label patchID) const;
			//*************** END based on solidGeneral*****************
			
            //- Return if this patch is the slave
            Switch master() const
            {
                return master_;
            }

            //- Return the shadow patch names
            //const wordList& shadowPatchNames() const;
			
			//*************** based on solidGeneral*****************
			//- Return the list of shadow patch names
            const List<word>& shadowPatchNames() const;
			//*************** END based on solidGeneral*****************

            //- Return the shadow patch indices
            const labelList& shadowPatchIndices() const;

            //- Return the shadow patch
            //  This function is not allowed if shadowPatchIndices.size() > 1
            const solid4GeneralContactFvPatchVectorField& shadowPatchField() const;

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

            //- Return const reference to deformed shadow zones
            const PtrList<globalPolyPatch>& shadowZones() const;

            //- Return reference to deformed shadow zones
            PtrList<globalPolyPatch>& shadowZones();

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

            //- Return reference to contactPerShadow fields
            PtrList<scalarField>& contactPerShadow();

            //- Return const reference to contactPerShadow fields
            const PtrList<scalarField>& contactPerShadow() const;


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
