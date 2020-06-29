/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    solidGenContactFvPatchVectorField

Description
    contact boundary conditions for finite volume solidMechanics solvers.
    There is a master and a slave patch. The master creates the contactModel
    and the slave gets the pointer to it.
    The treatment of the normal contact and tangential contact can be specified
    individually at run-time using the normalContactModel and
    frictionContactModel.

SourceFiles
    solidContactFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef solidContactFvPatchVectorField_H
#define solidContactFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "nonLinearGeometry.H"
#include "normalContactModel.H"
#include "frictionContactModel.H"
#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"
#include "ggiInterpolation.H"
#include "Switch.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidContactFvPatchVectorField
:
    public directionMixedFvPatchVectorField   //non-orthogonal (grid points) correction for diffusion term 
{
    // Private data
	
	//* * * * * * Added from the solidGeneralContactFvPatchVectorField.H  * * * * * * * * //
	
		//- Is this the global master patch: the global master patches keeps a
        // copy of mesh.allPoints() and updates them to the deformed position
        // each time before the contact is corrected
        mutable bool* globalMasterPtr_;
		
		//- Patch index of global master patch
        mutable label* globalMasterIndexPtr_;
		
		//- If the current patch is the local slave of the contact pair
        // The local slave is the patch with the lower patch index
        mutable boolList* localSlavePtr_;
		
		//- List of shadow patch names
        mutable List<word>* shadowPatchNamesPtr_;

        //- List of shadow patch indices
        mutable List<label>* shadowPatchIndicesPtr_;
		

        // Copy of contact dictionary
        const dictionary dict_;

        //- List of normal contact models
        mutable PtrList<normalContactModel> normalModels_;

        //- List of friction contact models
        mutable PtrList<frictionContactModel> frictionModels_;

        // Master face zone kept in deformed configuration
        mutable standAlonePatch* zonePtr_;

        //- List of zone-to-zone interpolators
        mutable PtrList<newGgiStandAlonePatchInterpolation> zoneToZones_;

        //- Global face zones in the mesh
        mutable labelList* globalFaceZonesPtr_;

        // Point map from global to local for global face zones
        mutable IOList<labelList>* globalToLocalFaceZonePointMapPtr_;

        //- The current patch tractions for each contatc pair
        mutable List<vectorField>* curPatchTractionPtr_;

        // Rate of friction energy dissipation
        mutable scalarField* QcPtr_;

        // Rate of friction energy dissipation for each contact pair
        mutable List<scalarField>* QcsPtr_;

        //- Offset used if bounding box of patch has a zero dimension
        mutable scalar bbOffset_;
		
	// * * * * * * * * * * * * * * * * * * * * * * End General * * * * * * * * * * * * * * * * * * * * //	

        //- Name of the displacement field
        const word fieldName_;

        //- is this patch the master
        const Switch master_;

        // is the contact boundary active
        Switch contactActive_;

        // is the contact boundary active //- Is the current patch considered rigid in contact
        Switch rigidMaster_;

        //- contact normallaw pointer
        normalContactModel* normalContactModelPtr_;
        frictionContactModel* frictionContactModelPtr_;

        //- shadow patch ID
        const label shadowPatchID_;

        //- either the mater of slave is in charge of correcting the
        // contact laws
        //- the patch with the lowest patchID is in charge of correcting
        // the contact laws
        const Switch patchInChargeOfCorrection_;

        // master and slave face zone names and IDs
        const word masterFaceZoneName_;
        const word slaveFaceZoneName_;
        const label masterFaceZoneID_;
        const label slaveFaceZoneID_;

        // the master will keep a copy of two global master and slave patches
        // these patches will be moved to the deformed position prior to
        // correction of the contact, and these patches are used to calculate
        // contact gaps, interpolation of tractions, etc.
        // we keep a copy of the points and faces and move them
        // pointField masterFaceZonePatchPoints_;
        // pointField slaveFaceZonePatchPoints_;
        // faceList masterFaceZonePatchFaces_;
        // faceList slaveFaceZonePatchFaces_;
        // PrimitivePatch<face, Foam::List, const pointField&>*
        // masterFaceZonePatchPtr_;
        // PrimitivePatch<face, Foam::List, const pointField&>*
        // slaveFaceZonePatchPtr_;
        // these patches keep a local copy of the points and faces
        PrimitivePatch<face, Foam::List, pointField>* masterFaceZonePatchPtr_;
        PrimitivePatch<face, Foam::List, pointField>* slaveFaceZonePatchPtr_;

        // method to interpolate traction from slave to master
        word interpolationMethod_;

       //- zoneToZone or ggiZone for interpolation of traction from slave to
       // master
       //zoneToZoneInterpolation* slaveToMasterPatchToPatchInterpolatorPtr_;
       PatchToPatchInterpolation< PrimitivePatch<face, Foam::List, pointField>,
                                  PrimitivePatch<face, Foam::List, pointField>
                                  >* slaveToMasterPatchToPatchInterpolatorPtr_;
       // we call it slaveToMasterGgi but it can interpolate either direction
       // ggiZoneInterpolation* slaveToMasterGgiInterpolatorPtr_;
       GGIInterpolation< PrimitivePatch< face, Foam::List, pointField >,
                         PrimitivePatch< face, Foam::List, pointField >
                         >* slaveToMasterGgiInterpolatorPtr_;

       // To interpolate value for the face zone patch faces to vertices
       // and vice versa
       PrimitivePatchInterpolation< PrimitivePatch<face, Foam::List, pointField>
                                    >* masterFaceZonePatchInterpolatorPtr_;
       PrimitivePatchInterpolation< PrimitivePatch<face, Foam::List, pointField>
                                    >* slaveFaceZonePatchInterpolatorPtr_;

        //- old face zone points
        vectorField oldMasterFaceZonePoints_;
        vectorField oldSlaveFaceZonePoints_;

        //- patch to patch algorithm
        const intersection::algorithm alg_;

        //- patch to patch direction
        const intersection::direction dir_;

       // current time index
       label curTimeIndex_;

       // current momentum iterator
       label iCorr_;

       // freqency the contact is corrected
       label correctionFreq_;

       // if the solver uses an orthotropic approach
       Switch orthotropic_;

       // field for visualisation of stick-slip regions
       const volScalarField* stickSlipFieldPtr_;

       // Force correction of contact on next call
       bool forceCorrection_;

      //- Is it a non linear solver
      nonLinearGeometry::nonLinearType nonLinear_;


    // Private Member Functions
		
		
	//* * * * * * Added from the solidGeneralContactFvPatchVectorField.H  * * * * * * * * //	
		
		//- Is a moving mesh (e.g. updated Lagrangian) approach be used
        bool movingMesh() const;

        //- Move the master and slave zones to the deformed configuration
        void moveZonesToDeformedConfiguration();

        //- Check global master patch
        void calcGlobalMaster() const;

        //- Check global master patch index
        void calcGlobalMasterIndex() const;

        //- Calculate local slaves
        void calcLocalSlave() const;

        //- Return local slaves
        const boolList& localSlave() const;

        //- Calculate shadow patch names
        void calcShadowPatchNames() const;

        //- Calculate shadow zone names
        void calcShadowZoneNames() const;

        //- Create the normal contact model
        void calcNormalModels() const;

        //- Create the friction contact model
        void calcFrictionModels() const;

        //- Calculate zone index
        void calcZoneIndex() const;

        //- Calculate zone
        void calcZone() const;

        //- Calculate zone old points
        void calcOldZonePoints() const;

        //- Calculate interpolation
        void calcZoneToZones() const;

        //- Make the current patch traction list
        void makeCurPatchTractions() const;

        //- Calculate rate of dissipated friction energy
        void calcQc() const;

        //- Calculate rate of dissipated friction energy for each contact
        //  pair
        void calcQcs() const;

        //- Calculate the bounding box offset
        void calcBbOffset() const;

        //- Return the bounding box offset
        scalar bbOffset() const;

        // Contact dictionary
        const dictionary& dict() const
        {
            return dict_;
        };

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
	
	// * * * * * * * * * * * * * * * * * * * * * * End General * * * * * * * * * * * * * * * * * * * * //
	
		//  Move the contact face zone patches to the deformed position
        virtual void moveFaceZonePatches();

        //  check shadow patch and face zones exist
        virtual bool checkPatchAndFaceZones(const dictionary& dict) const;

public:

	//* * * * * * Added from the solidGeneralContactFvPatchVectorField.H  * * * * * * * * //
	
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

        //- Construct by mapping given solidGeneralContactFvPatchVectorField
        // onto a new patch
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
	
	// Destructor

        virtual ~solidGeneralContactFvPatchVectorField();
	
	
	// * * * * * * * * * * * * * * * * * * * * * * End General * * * * * * * * * * * * * * * * * * * * //



    // Destructor

//         virtual ~solidContactFvPatchVectorField();


    // Member functions
	
		//* * * * * * Added from the solidGeneralContactFvPatchVectorField.H  * * * * * * * * //
		
		
		// Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }

            //- Return true if this patch is the global master
            bool globalMaster() const;

            //- Return index of global master patch
            label globalMasterIndex() const;

            //- Return the list of shadow patch names
            const List<word>& shadowPatchNames() const;

            //- Return the list of shadow patch indices
            const List<label>& shadowPatchIndices() const;

            //- Return index of the zone
            label zoneIndex() const;

            //- Return the list of shadow zone names
            const List<word>& shadowZoneNames() const;

            //- Return the list of shadow zone indices
            const List<label>& shadowZoneIndices() const;

            //- Return reference to the normal contact model
            normalContactModel& normalModel(const label shadowI);

            //- Return const reference to the normal contact model
            const normalContactModel& normalModel(const label shadowI) const;

            //- Return reference to the friction contact model
            frictionContactModel& frictionModel(const label shadowI);

            //- Return const reference to the friction contact model
            const frictionContactModel& frictionModel
            (
                const label shadowI
            ) const;

            //- Return const reference to deformed master zone
            const standAlonePatch& zone() const;

            //- Return reference to deformed master zone
            standAlonePatch& zone();

            //- Return const reference to deformed master zone
            const standAlonePatch& shadowZone(const label shadowI) const;

            //- Return reference to deformed master zone
            standAlonePatch& shadowZone(const label shadowI);

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

            //- Given the patch index of a shadow, return the position in the
            // shadow list of the current patch
            label findShadowID(const label patchID) const;

            //- Return the current patch traction for shadowI
            const vectorField& curPatchTractions(const label shadowI) const;

            //- Return the current patch traction for shadowI
            vectorField& curPatchTractions(const label shadowI);
	
		// * * * * * * * * * * * * * * * * * * * * * * End General * * * * * * * * * * * * * * * * * * * * //

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }

            //- Return if this patch is the master
            Switch master() const
            {
                return master_;
            }

            //- Return contact normal law pointer
            normalContactModel* normalContactModelPtr() const
            {
                return normalContactModelPtr_;
            }

            //- Return contact normal friction pointer
            frictionContactModel* frictionContactModelPtr() const
            {
                return frictionContactModelPtr_;
            }

            //- Return masterFaceZonePatch pointer
            PrimitivePatch<face, Foam::List, pointField>*
            masterFaceZonePatchPtr() const
            {
                return masterFaceZonePatchPtr_;
            }

            //- Return masterFaceZonePatch pointer
            PrimitivePatch<face, Foam::List, pointField>*
            slaveFaceZonePatchPtr() const
            {
                return slaveFaceZonePatchPtr_;
            }

            //- Return patchToPatch interpolator pointer
            PatchToPatchInterpolation
            <
                PrimitivePatch<face, Foam::List, pointField>,
                PrimitivePatch<face, Foam::List, pointField>
                >* slaveToMasterPatchToPatchInterpolatorPtr() const
            {
                return slaveToMasterPatchToPatchInterpolatorPtr_;
            }

            //- Return GGI interpolator pointer
            GGIInterpolation
            <
                PrimitivePatch< face, Foam::List, pointField >,
                PrimitivePatch< face, Foam::List, pointField >
                >* slaveToMasterGgiInterpolatorPtr() const
                {
                    return slaveToMasterGgiInterpolatorPtr_;
                }

            //- Return correction frequency
            label correctionFreq() const
            {
                return correctionFreq_;
            }

            //  Return contactActive
            virtual Switch contactActive() const
            {
                return contactActive_;
            }

        // Mapping functions
		
		//* * * * * * Added from the solidGeneralContactFvPatchVectorField.H  * * * * * * * * //
		
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
		
		
		// * * * * * * * * * * * * * * * * * * * * * * End General * * * * * * * * * * * * * * * * * * * * //

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
		
		//* * * * * * Added from the solidGeneralContactFvPatchVectorField.H  * * * * * * * * * *//
		
		    //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
		
		// * * * * * * * * * * * * * * * * * * * * * * End General * * * * * * * * * * * * * * * * * * * * //

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
                );

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;


        // Member functions
		
		//* * * * * * Added from the solidGeneralContactFvPatchVectorField.H  * * * * * * * * * *//
		
		    //- Dissipated energy due to friction
            virtual const scalarField& Qc() const;

            //- Dissipated energy due to friction from contact with shadowI
            virtual const scalarField& Qc(const label shadowI) const;

            //- Write
            virtual void write(Ostream&) const;
		
		// * * * * * * * * * * * * * * * * * * * * * * End General * * * * * * * * * * * * * * * * * * * * //

            // interpolate the given field from patch to patch in parallel
            virtual tmp<vectorField> interpolateSlaveToMaster
            (
                const vectorField slaveField
                );

            //- Force contact correction
            virtual void forceCorrection()
            {
                forceCorrection_ = true;
            }

            //- Dissipated energy due to friction
            virtual tmp<scalarField> Qc() const;

            // Return const reference to stickSlipField
            const scalarField& stickSlipField() const
            {
                if (!stickSlipFieldPtr_)
                {
                    FatalError << "solidContact::stickSlipField() pointer"
                               << " is NULL" << exit(FatalError);
                }

                return stickSlipFieldPtr_->boundaryField()[patch().index()];
            }

            //- Write
            virtual void write(Ostream&) const;

    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
