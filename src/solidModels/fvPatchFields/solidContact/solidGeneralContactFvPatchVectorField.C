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

\*---------------------------------------------------------------------------*/

#include "solidGeneralContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// *************************************** START general ****************************

bool Foam::solidGeneralContactFvPatchVectorField::movingMesh() const
{
    // If the deformation gradient "F" and the displacement increment DU" are
    // found then we can assume it is a moving mesh (updated Lagrangian) case
    if
    (
        db().foundObject<volVectorField>("DU")
     && db().foundObject<volTensorField>("F")
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::solidGeneralContactFvPatchVectorField::
moveFaceZonesToDeformedConfiguration()   // CHECK ONLY in Deformed Configuration 
{
	Info<<"Step1: In moveFaceZonesToDeformedConfiguration():"<<__LINE__<<endl;
    Info<<"The contact surfaces meshes are moved to their deformed position "<<__LINE__<<endl;
    // Only the master moves the zones
    if (!globalMaster())
    {
        return;
    }

    // Reference to mesh for tidiness
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Method
    // We will interpolate the patch face displacements to the patch vertices
    // and then add these vertex/point displacements to the initial patch
    // points
    // We need to take care in parallel, and also realise that the solidModel
    // might have a moving or stationary mesh
	
    // Shadow patch and zone indices
    const labelList& shadPatchIndices = shadowPatchIndices();
    const labelList& shadZoneIndices = shadowZoneIndices();
	

    forAll(shadPatchIndices, shadowI)
    {
        // Assemble the zone face displacement field to move the zones
        vectorField zoneD(zone().size(), vector::zero);
        vectorField shadowZoneD(shadowZone(shadowI).size(), vector::zero);

        // For a non-moving mesh, we will move the zones by the total
        // displacement, whereas for a moving mesh (updated Lagrangian), we will
        // move the zones by the displacement increment

        if (movingMesh())
        {
            // Updated Lagrangian, so we will move the zones by the displacement
            // increment

            // Lookup the current total displacement field
            const volVectorField& DD = db().lookupObject<volVectorField>("DU");

            // Take a reference to the patch face displacement increment field
            const vectorField& patchDD =
                DD.boundaryField()[patch().index()];
            const vectorField& shadowPatchDD =
                DD.boundaryField()[shadPatchIndices[shadowI]];

            zoneD =
                zoneField(zoneIndex(), patch().index(), patchDD);
            shadowZoneD =
                    zoneField
                    (
                        shadZoneIndices[shadowI],
                        shadPatchIndices[shadowI],
                        shadowPatchDD
                    );
        }
        else
        {
            // Non-moving mesh: we will move the zones by the total displacement

            // Lookup the current total displacement field
            const volVectorField& D = db().lookupObject<volVectorField>("U");

            // Take a reference to the patch face total displacement field
            const vectorField& patchD =
                D.boundaryField()[patch().index()];

            //const vectorField& shadowPatchDD = 
			const vectorField& shadowPatchD = 
                D.boundaryField()[shadPatchIndices[shadowI]];

            zoneD =
                zoneField(zoneIndex(), patch().index(), patchD);
            shadowZoneD =
                zoneField
                (
                    shadZoneIndices[shadowI],
                    shadPatchIndices[shadowI],
                    shadowPatchD   //shadowPatchDD 
                );
        }

        // Interpolate the zone face field to the zone points
        const pointField zonePointD =
            zoneFaceToPointInterpolate(zoneIndex(), zoneD, -1);
        const pointField shadowZonePointD =
            zoneFaceToPointInterpolate
            (
                shadZoneIndices[shadowI], shadowZoneD, shadowI
            );
		
		// The zone deformed points are the initial position plus the
        // displacement
        const pointField zoneNewPoints =
            mesh.faceZones()[zoneIndex()]().localPoints()
          + zonePointD;
        const pointField shadowZoneNewPoints =
            mesh.faceZones()[shadZoneIndices[shadowI]]().localPoints()
          + shadowZonePointD;

        // Move the zones

        // Remove zones weights
        if (shadowI == 0)
        {
			Info<<"In moveFaceZonesToDeformedConfiguration():"<<__LINE__<<endl;
            zone().movePoints(zoneNewPoints);
        }
        shadowZone(shadowI).movePoints(shadowZoneNewPoints);

        // We need to use const_cast to move the standAlonePatch points as the
        // movePoints function only clears weights
        // Also, be careful to move the points are opposed to the localPoints
        if (shadowI == 0)
        {
            const_cast<pointField&>(zone().points()) = zoneNewPoints;
        }
        const_cast<pointField&>(shadowZone(shadowI).points()) =
            shadowZoneNewPoints;
    }
}

void Foam::solidGeneralContactFvPatchVectorField::calcGlobalMaster() const   //// CHECK method to calculate global master
{
    if (globalMasterPtr_)
    {
        FatalErrorIn
            (
                "void Foam::solidGeneralContactFvPatchVectorField::"
                "calcGlobalMaster() const"
            )   << "globalMasterPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solidGeneralContact patch i.e. the one
    // with the lowest patch index
	
	Info<< "In calcGlobalMaster() "<<__LINE__<<endl;
	Info<< "globalMasterIndex() in calcGlobalMaster() "<<globalMasterIndex()<<endl;
	Info<< "patch().index() in calcGlobalMaster() "<<patch().index()<<endl;
    Info<< "patch().name() in calcGlobalMaster() "<<patch().name()<<endl;
	if (globalMasterIndex() == patch().index())
    {
        globalMasterPtr_ = new bool(true);
    }
    else
    {
        globalMasterPtr_ = new bool(false);
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcGlobalMasterIndex() const
{
    if (globalMasterIndexPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcGlobalMasterIndex() const"
        )   << "globalMasterIndexPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solidGeneralContact patch i.e. the one
    // with the lowest patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    globalMasterIndexPtr_ = new label(-1);
    label& gMasterID = *globalMasterIndexPtr_;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solidGeneralContactFvPatchVectorField::typeName
        )
        {
            gMasterID = patchI;

            break;
        }
    }

    // Check there is only one global master

    label GMasterID = returnReduce(gMasterID, maxOp<label>());

    if (gMasterID != GMasterID)
    {
        FatalErrorIn
        (
            "solidGeneralContactFvPatchVectorField::"
            "calcGlobalMasterIndex() const"
        )   << "There are multiple global masters" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< nl << "The global master contact patch is "
            << patch().boundaryMesh()[gMasterID].name() << endl;
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcLocalSlave() const
{
    if (localSlavePtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcLocalSlave() const"
        )   << "localSlavePtr_ already set" << abort(FatalError);
    }

    localSlavePtr_ = new boolList(shadowPatchNames().size(), false);

    boolList& localSlave = *localSlavePtr_;
	
	Info<< "In calcLocalSlave() "<<__LINE__<<endl;
	Info<< "patch().index() in calcLocalSlave() "<<patch().index()<<endl;
    forAll(localSlave, shadowI)
    {
        if (patch().index() < shadowPatchIndices()[shadowI])
        {
            localSlave[shadowI] = true;

            Info<< "solidGeneralContact: "
                << shadowPatchNames()[shadowI] << " (master)" << " to "
                << patch().name() << " (slave)" << endl;
        }
    }
}

const Foam::boolList&
Foam::solidGeneralContactFvPatchVectorField::localSlave() const
{
    if (!localSlavePtr_)
    {
        calcLocalSlave();
    }
	Info<< "The current field in localSlave() is "<< dimensionedInternalField().name()<< endl;
	Info<< "The current patch in localSlave() is "<< patch().name()<< endl;			
	Info<< "The size of current patch in localSlave() is "<< patch().size()<< endl;
	Info<< "*localSlavePtr_ in localSlave() is "<<*localSlavePtr_<< endl;
    return *localSlavePtr_;
}

void Foam::solidGeneralContactFvPatchVectorField::calcShadowPatchNames() const
{
    if (shadowPatchNamesPtr_ || shadowPatchIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcShadowPatchNames() const"
        )   << "shadowPatchNames_ or shadowPatchIndices_ already set"
            << abort(FatalError);
    }

    // Add each solidGeneralContact patch in the order of increasing patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    // Count shadow patches

    label nShadPatches = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solidGeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
            nShadPatches++;
        }
    }
	
    shadowPatchNamesPtr_ = new wordList(nShadPatches);
    wordList& shadowPatchNames = *shadowPatchNamesPtr_;
	Info<<"In calcShadowPatchNames():"<<__LINE__<<endl;
//	Info<<"What is *shadowPatchNamesPtr_? in calcShadowPatchNames():"<<*shadowPatchNamesPtr_<<endl;

    shadowPatchIndicesPtr_ = new labelList(nShadPatches);
    labelList& shadowPatchIndices = *shadowPatchIndicesPtr_;
//	Info<<"In calcShadowPatchNames():"<<__LINE__<<endl;
	
    // Record shadow patch names

    label shadowI = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solidGeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
            shadowPatchNames[shadowI] = patch().boundaryMesh()[patchI].name();

            shadowPatchIndices[shadowI++] = patchI;
        }
    }
}



 void Foam::solidGeneralContactFvPatchVectorField::calcShadowZoneNames() const
{
    if (shadowZoneNamesPtr_ || shadowZoneIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcShadowZoneNames() const"
        )   << "shadowZoneNames_ or shadowZoneIndices_ already set"
            << abort(FatalError);
    }

    const wordList& shadNames = shadowPatchNames();

    shadowZoneNamesPtr_ = new wordList(shadNames.size());
    wordList& shadowZoneNames = *shadowZoneNamesPtr_;

    shadowZoneIndicesPtr_ = new labelList(shadNames.size());
    labelList& shadowZoneIndices = *shadowZoneIndicesPtr_;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    forAll(shadNames, shadowI)
    {
        word zoneName = shadNames[shadowI] + "FaceZone";

        faceZoneID zone(zoneName, mesh.faceZones());

        if (!zone.active())
        {
            FatalErrorIn("solidGeneralContactFvPatchVectorField")
                << "Face zone name " << zoneName
                << " not found.  Please check your zone definition."
                << abort(FatalError);
        }

        shadowZoneNames[shadowI] = zoneName;

        shadowZoneIndices[shadowI] = zone.index();
    }
}

void Foam::solidGeneralContactFvPatchVectorField::calcNormalModels() const
{
    if (!normalModels_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "calcNormalModel() const"
        )   << "normalModels already set" << abort(FatalError);
    }

    normalModels_.setSize(shadowPatchNames().size());

    const boolList& locSlave = localSlave();

    forAll(normalModels_, shadowI)
    {
        // Only the local slave creates the contact model
        if (locSlave[shadowI])
        {
			
            // Calculate normal contact forces
            normalModels_.set
            (
                shadowI,
                generalNormalContactModel::New
                (
                    word(dict().lookup("generalNormalContactModel")),
                    patch().boundaryMesh()[shadowPatchIndices()[shadowI]],
                    dict(),
                    shadowPatchIndices()[shadowI], // master
                    patch().index(), // slave
                    shadowZone(shadowI), // master
                    zone() // slave
                )
            );
			
        }
    }

}

Foam::generalNormalContactModel&
Foam::solidGeneralContactFvPatchVectorField::normalModel(const label shadowI)
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("normalModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (normalModels_.empty())
    {
        calcNormalModels();
    }

    return normalModels_[shadowI];
}

const Foam::generalNormalContactModel&
Foam::solidGeneralContactFvPatchVectorField::normalModel
(
    const label shadowI
) const
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("normalModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (normalModels_.empty())
    {
        calcNormalModels();
    }

    return normalModels_[shadowI];
}


void Foam::solidGeneralContactFvPatchVectorField::calcFrictionModels() const
{
  
  if (!frictionModels_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "calcFrictionModel(shadowI) const"
        )   << "frictionModelPtr_[shadowI] already set" << abort(FatalError);
    }

    frictionModels_.setSize(shadowPatchNames().size());

    const boolList& locSlave = localSlave();

    forAll(frictionModels_, shadowI)
    {
        if (locSlave[shadowI])
        {
            frictionModels_.set
                (
                    shadowI,
                    generalFrictionContactModel::New
                    (
                        word(dict().lookup("generalFrictionContactModel")),
                        patch().boundaryMesh()[shadowPatchIndices()[shadowI]],
                        dict(),
                        shadowPatchIndices()[shadowI], // master
                        patch().index() // slave
                    )
                );
        }
    }
	
}

Foam::generalFrictionContactModel&
Foam::solidGeneralContactFvPatchVectorField::frictionModel(const label shadowI)
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("frictionModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (frictionModels_.empty())
    {
        calcFrictionModels();
    }

    return frictionModels_[shadowI];
}


const Foam::generalFrictionContactModel&
Foam::solidGeneralContactFvPatchVectorField::frictionModel
(
    const label shadowI
) const
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("frictionModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (frictionModels_.empty())
    {
        calcFrictionModels();
    }

    return frictionModels_[shadowI];
}

void Foam::solidGeneralContactFvPatchVectorField::calcZoneIndex() const
{
    if (zoneIndex_ != -1)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcZoneIndex()"
            "const"
        )   << "zoneIndex_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    word zoneName = patch().name() + "FaceZone";

    faceZoneID zone(zoneName, mesh.faceZones());

    if (!zone.active())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcZoneIndex()"
            "const"
        )   << "Face zone name " << zoneName
            << " not found.  Please check your zone definition." << nl
            << "Current faceZones are:" << mesh.faceZones().names()
            << abort(FatalError);
    }

    zoneIndex_ = zone.index();
}

void Foam::solidGeneralContactFvPatchVectorField::calcZone() const
{
    if (zonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcZone() const"
        )   << "zonePtr_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    zonePtr_ =
        new standAlonePatch
        (
            mesh.faceZones()[zoneIndex()]().localFaces(),
            mesh.faceZones()[zoneIndex()]().localPoints()
        );
}

//**************************start definition from solid4Foam **********************************************

void Foam::solidGeneralContactFvPatchVectorField::calcShadowZonesNewGgi() const
{
	Info<<"CHECK solid4Foam Here I am in calcShadowZonesNewGgi()"<<__LINE__<<endl;
    if (debug)
    {
        InfoIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcShadowZonesNewGgi() const"
        )   << patch().name() << " : making the shadow zones" << endl;
    }
	
	const boolList& locSlave = localSlave();

//************************ start ERROR - skip this shadowI not defined ***********

//    if (!locSlave[shadowI])
	if (!globalMasterPtr_)
//	if (!master_)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "calcShadowZonesNewGgi() const"
        )   << "Trying to create shadow zones on a slave" << abort(FatalError);
    }
	
//************************ end ERROR - skip this shadowI not defined ***********

    if (!shadowZonesNewGgi_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "calcShadowZonesNewGgi() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const wordList& shadPatchNames = shadowPatchNames();

    shadowZonesNewGgi_.setSize(shadPatchNames.size());

    forAll(shadowZonesNewGgi_, shadPatchI)
    {
		Info<<"Here I am in calcShadowZonesNewGgi()"<<__LINE__<<endl;
		Info<<"shadPatchNames[shadPatchI] in calcShadowZonesNewGgi()"<<shadPatchNames[shadPatchI]<<endl;
        // Note: the main mesh will either be in the initial configuration or
        // the updated configuration
        shadowZonesNewGgi_.set
        (
            shadPatchI,
            new globalPolyPatch
            (
                shadPatchNames[shadPatchI],
                patch().boundaryMesh().mesh()
            )
        );
    }
	
}

//**************************END definition from solid4Foam **********************************************

Foam::label Foam::solidGeneralContactFvPatchVectorField::zoneIndex() const
{
	Info<<"Does it enter here? in zoneIndex()"<<__LINE__<<endl;	
    if (zoneIndex_ == -1)
    {
        calcZoneIndex();
    }
	Info<<"Does it enter here? in zoneIndex()"<<__LINE__<<endl;
	Info<<"zoneIndex_ in zoneIndex()"<<zoneIndex_<<endl;
    return zoneIndex_;
}

const Foam::standAlonePatch&
Foam::solidGeneralContactFvPatchVectorField::zone() const
{
    if (!zonePtr_)
    {
        calcZone();
    }
	Info<< "The current field in zone() is "<< dimensionedInternalField().name()<< endl;
//	Info<< "The current dimensionedInternalField().size() in zone() is "<< dimensionedInternalField().size()<< endl;
	Info<< "The current patch in zone() is "<< patch().name()<< endl;			
	Info<< "The size of current patch in zone() is "<< patch().size()<< endl;
	Info<< "*zonePtr_ in zone() is "<<*zonePtr_<< endl;
    return *zonePtr_;
}


Foam::standAlonePatch& Foam::solidGeneralContactFvPatchVectorField::zone()
{
    if (!zonePtr_)
    {
        calcZone();
    }

    return *zonePtr_;
}

const Foam::standAlonePatch&
Foam::solidGeneralContactFvPatchVectorField::shadowZone
(
    const label shadowI
) const
{
    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    const solidGeneralContactFvPatchVectorField& shadowPatchField =
        refCast<const solidGeneralContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndices()[shadowI]]
        );

    return shadowPatchField.zone();
}


Foam::standAlonePatch&
Foam::solidGeneralContactFvPatchVectorField::shadowZone
(
    const label shadowI
)
{
    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    // Const cast away the const-ness
    solidGeneralContactFvPatchVectorField& shadowPatchField =
        const_cast<solidGeneralContactFvPatchVectorField&>
        (
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadowI]]
            )
        );

    return shadowPatchField.zone();
}

//Here I am update 


void Foam::solidGeneralContactFvPatchVectorField::calcQc() const
{
    if (QcPtr_)
    {
        FatalErrorIn("solidGeneralContactFvPatchVectorField::calcQc")
            << "QcPtr_ already set!" << abort(FatalError);
    }

    QcPtr_ = new scalarField(patch().size(), 0.0);

    scalarField& Qc = *QcPtr_;

    // For now, we assume traction is constant over time-step
    // Todo: use trapezoidal rule
    vectorField curTraction(Qc.size(), vector::zero);

    // sigma/sigmaCauchy is up-to-date as Qc is called after momentum loop
    // has converged and sigma has been updated and mesh moved
    if
    (
        db().objectRegistry::foundObject<volSymmTensorField>
        (
            "sigmaCauchy"
        )
    )
    {
        const symmTensorField& sigma =
            db().objectRegistry::lookupObject<volSymmTensorField>
            (
                "sigmaCauchy"
            ).boundaryField()[patch().index()];

        curTraction = patch().nf() & sigma;
    }
    else
    {
        const symmTensorField& sigma =
            db().objectRegistry::lookupObject<volSymmTensorField>
            (
                "sigma"
            ).boundaryField()[patch().index()];

        curTraction = patch().nf() & sigma;
    }


    // Accumulate Qc for from all shadows

    const scalar deltaT = patch().boundaryMesh().mesh().time().deltaT().value();

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    const boolList& locSlave = localSlave();

    forAll(locSlave, shadowI)
    {
        vectorField curPatchSlip(patch().size(), vector::zero);

        // Calculate slip
        if (locSlave[shadowI])
        {
			
            curPatchSlip = frictionModel(shadowI).slip();
			
        }
        else
        {
            const solidGeneralContactFvPatchVectorField& shadowPatchField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadowI]]
            );

            const label locShadowID =
                shadowPatchField.findShadowID(patch().index());
			
			
            
			vectorField shadowPatchSlip =
                shadowPatchField.frictionModel(locShadowID).slip();
				
			
			
            
			vectorField shadowZoneSlip =
                zoneField
                (
                    shadowZoneIndices()[shadowI],
                    shadowPatchIndices()[shadowI],
                    shadowPatchSlip
                );
			

            // Interpolate from shadow to the current patch
            // Face-to-face
			
			
			vectorField curZoneSlip =
                shadowPatchField.zoneToZoneNewGgi(locShadowID).slaveToMaster
                (
                    shadowZoneSlip
                );
			
			
            
			curPatchSlip =
                patchField
                (
                    patch().index(),
                    zoneIndex(),
                    curZoneSlip
                );
				
			
        }

        // Heat flux rate: rate of dissipated frictional energy
        // The dot product of the traction vectors and the slip vectors
        // gives the dissipated frictional energy rate per unit area, which
        // is always positive
        Qc += mag(curTraction & (curPatchSlip/deltaT));
    }
}



 
// *************************************** END general ****************************

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidGeneralContactFvPatchVectorField::solidGeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:	// *************************************** START general ****************************
	solidTractionFvPatchVectorField(p, iF),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
	dict_(NULL),
	normalModels_(0),
    frictionModels_(0),
	zonePtr_(NULL),
	zoneIndex_(-1),	
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL), 
	zoneToZones_(0),
	zoneToZonesNewGgi_(0),
	rigidMaster_(false),
	curPatchTractionPtr_(NULL),
	alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
	curTimeIndex_(-1),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(0.0)
	
	// *************************************** END general ****************************
//    shadowPatchID_(-1),

//    masterFaceZoneID_(-1),
//    slaveFaceZoneID_(-1)

// ******************************************** START General *****************************************
{
	Info<<"Does it enter here? in C1(p, iF)"<<__LINE__<<endl;
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        ) << endl;
    }
}
// ********************************************** END General ********************************************

solidGeneralContactFvPatchVectorField::solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:	// ******************************************** START General *****************************************
	
	solidTractionFvPatchVectorField(ptf, p, iF, mapper),
	dict_(ptf.dict_),
	globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),  // (NULL),
    frictionModels_(ptf.frictionModels_), // (NULL),
    zonePtr_(NULL),
    zoneToZones_(0),
	zoneToZonesNewGgi_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)
// ********************************************** END General ********************************************


 //   shadowPatchID_(ptf.shadowPatchID_),
	
 //   masterFaceZoneID_(ptf.masterFaceZoneID_),
 //   slaveFaceZoneID_(ptf.slaveFaceZoneID_)
	
// ******************************************** START General *****************************************
{
	Info<<"Does it enter here? in C2(ptf, p, iF, mapper)"<<__LINE__<<endl;
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF,"
            "    const fvPatchFieldMapper& mapper"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.globalMasterPtr_)
    {
        globalMasterPtr_ = new bool(*ptf.globalMasterPtr_);
    }

    if (ptf.globalMasterIndexPtr_)
    {
        globalMasterIndexPtr_ = new label(*ptf.globalMasterIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.shadowPatchNamesPtr_)
    {
        shadowPatchNamesPtr_ = new wordList(*ptf.shadowPatchNamesPtr_);
    }

    if (ptf.shadowPatchIndicesPtr_)
    {
        shadowPatchIndicesPtr_ = new labelList(*ptf.shadowPatchIndicesPtr_);
    }

    if (ptf.shadowZoneNamesPtr_)
    {
        shadowZoneNamesPtr_ = new wordList(*ptf.shadowZoneNamesPtr_);
    } 

    if (ptf.shadowZoneIndicesPtr_)
    {
        shadowZoneIndicesPtr_ = new labelList(*ptf.shadowZoneIndicesPtr_);
    }	

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (!ptf.zoneToZonesNewGgi_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solidGeneralContact: zoneToZone GGI interpolators not mapped"
            << endl;
    }

    if (ptf.curPatchTractionPtr_)
    {
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }

    if (ptf.QcPtr_)
    {
        QcPtr_ = new scalarField(*ptf.QcPtr_);
    }

    if (ptf.QcsPtr_)
    {
        QcsPtr_ = new List<scalarField>(*ptf.QcsPtr_);
    }
}

// ********************************************** END General ********************************************


solidGeneralContactFvPatchVectorField::solidGeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
    : 
	// ********************************************** START General ********************************************
	
	solidTractionFvPatchVectorField(p, iF),
	dict_(dict),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(-1),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(dict.lookupOrDefault<Switch>("rigidMaster", false)),
    normalModels_(0),
    frictionModels_(0),
    zonePtr_(0),
    zoneToZones_(0),
	zoneToZonesNewGgi_(0),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(0.0)
	
	// ********************************************** END General ********************************************

  /*  shadowPatchID_
    (
        patch().patch().boundaryMesh().findPatchID(dict.lookup("shadowPatch"))
        ),
	
    masterFaceZoneID_
    (
        patch().boundaryMesh().mesh().faceZones().findZoneID
        (
            masterFaceZoneName_
            )
        ),
    slaveFaceZoneID_
    (
        patch().boundaryMesh().mesh().faceZones().findZoneID(slaveFaceZoneName_)
        )
	*/
	
	// ********************************************** START General ********************************************
	
	{
		Info<<"Does it enter here? in C3(p, iF, dict)"<<__LINE__<<endl;
    Info<< "Creating " << solidGeneralContactFvPatchVectorField::typeName
        << " patch" << endl;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    } 
}
	
	// ********************************************** END General ********************************************


//******************************************** START General ******************************************** 
solidGeneralContactFvPatchVectorField::solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf
)
:
	solidTractionFvPatchVectorField(ptf),
	dict_(ptf.dict_),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
//	masterFaceZoneID_(0),
//	slaveFaceZoneID_(0),
    localSlavePtr_(NULL),
//	shadowPatchID_(0),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),   //(NULL),
    frictionModels_(ptf.frictionModels_),   //(NULL),
    zonePtr_(NULL),
    zoneToZones_(0),
	zoneToZonesNewGgi_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)
{
	Info<<"Does it enter here? in C4(ptf)"<<__LINE__<<endl;
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.globalMasterPtr_)
    {
        globalMasterPtr_ = new bool(*ptf.globalMasterPtr_);
    }

    if (ptf.globalMasterIndexPtr_)
    {
        globalMasterIndexPtr_ = new label(*ptf.globalMasterIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.shadowPatchNamesPtr_)
    {
        shadowPatchNamesPtr_ = new wordList(*ptf.shadowPatchNamesPtr_);
    }

    if (ptf.shadowPatchIndicesPtr_)
    {
        shadowPatchIndicesPtr_ = new labelList(*ptf.shadowPatchIndicesPtr_);
    }

    if (ptf.shadowZoneNamesPtr_)
    {
        shadowZoneNamesPtr_ = new wordList(*ptf.shadowZoneNamesPtr_);
    } 

    if (ptf.shadowZoneIndicesPtr_)
    {
        shadowZoneIndicesPtr_ = new labelList(*ptf.shadowZoneIndicesPtr_);
    } 

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (!ptf.zoneToZonesNewGgi_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solidGeneralContact: zoneToZone GGI interpolators not mapped"
            << endl;
    }

    if (ptf.curPatchTractionPtr_)
    {
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }

    if (ptf.QcPtr_)
    {
        QcPtr_ = new scalarField(*ptf.QcPtr_);
    }

    if (ptf.QcsPtr_)
    {
        QcsPtr_ = new List<scalarField>(*ptf.QcsPtr_);
    }
}

//**************************************************** END General**********************************************

solidGeneralContactFvPatchVectorField::solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
: 	//**************************************************** START General**********************************************

	solidTractionFvPatchVectorField(ptf, iF),
	dict_(ptf.dict_),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),  //(NULL),
    frictionModels_(ptf.frictionModels_), //(NULL),
    zonePtr_(NULL),
    zoneToZones_(0),
	zoneToZonesNewGgi_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
	curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)

	//**************************************************** END General**********************************************

/*	shadowPatchID_(ptf.shadowPatchID_),

    masterFaceZoneID_(ptf.masterFaceZoneID_),
    slaveFaceZoneID_(ptf.slaveFaceZoneID_) 
	*/

//**************************************************** START General**********************************************
{
	Info<<"Does it enter here? in C5(ptf, iF)"<<__LINE__<<endl;
	Info<<"In constructor 5 line :"<<__LINE__<<endl;
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.globalMasterPtr_)
    {
		Info<<"Does it enter here? in C5(ptf, iF)"<<__LINE__<<endl;
        globalMasterPtr_ = new bool(*ptf.globalMasterPtr_);
    }

    if (ptf.globalMasterIndexPtr_)
    {
		Info<<"Does it enter here? in C5(ptf, iF)"<<__LINE__<<endl;
        globalMasterIndexPtr_ = new label(*ptf.globalMasterIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
		Info<<"Does it enter here? in C5(ptf, iF)"<<__LINE__<<endl;
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.shadowPatchNamesPtr_)
    {
        shadowPatchNamesPtr_ = new wordList(*ptf.shadowPatchNamesPtr_);
    }

    if (ptf.shadowPatchIndicesPtr_)
    {
        shadowPatchIndicesPtr_ = new labelList(*ptf.shadowPatchIndicesPtr_);
    }

    if (ptf.shadowZoneNamesPtr_)
    {
        shadowZoneNamesPtr_ = new wordList(*ptf.shadowZoneNamesPtr_);
    } 

    if (ptf.shadowZoneIndicesPtr_)
    {
        shadowZoneIndicesPtr_ = new labelList(*ptf.shadowZoneIndicesPtr_);
    }	

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (!ptf.zoneToZonesNewGgi_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solidGeneralContact: zoneToZone GGI interpolators not mapped"
            << endl;
    }

    if (ptf.curPatchTractionPtr_)
    {
		Info<<"Does it enter here? in C5(ptf, iF)"<<__LINE__<<endl;
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }

    if (ptf.QcPtr_)
    {
        QcPtr_ = new scalarField(*ptf.QcPtr_);
    }

    if (ptf.QcsPtr_)
    {
        QcsPtr_ = new List<scalarField>(*ptf.QcsPtr_);
    }
}

//**************************************************** END General**********************************************

// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //

//**************************************************** START General**********************************************

Foam::solidGeneralContactFvPatchVectorField::
~solidGeneralContactFvPatchVectorField()
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "~solidGeneralContactFvPatchVectorField()"
        ) << endl;
    }

    deleteDemandDrivenData(globalMasterPtr_);
    deleteDemandDrivenData(globalMasterIndexPtr_);
    deleteDemandDrivenData(localSlavePtr_);
    deleteDemandDrivenData(shadowPatchNamesPtr_);
    deleteDemandDrivenData(shadowPatchIndicesPtr_);
    deleteDemandDrivenData(shadowZoneNamesPtr_);
    deleteDemandDrivenData(shadowZoneIndicesPtr_);

    normalModels_.clear();
    frictionModels_.clear();

    deleteDemandDrivenData(zonePtr_);

    zoneToZones_.clear();
	zoneToZonesNewGgi_.clear();

    deleteDemandDrivenData(curPatchTractionPtr_);
    deleteDemandDrivenData(QcPtr_);
    deleteDemandDrivenData(QcsPtr_);
}

//**************************************************** END General**********************************************

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

//**************************************************** START General**********************************************


Foam::label Foam::solidGeneralContactFvPatchVectorField::globalMasterIndex()
const
{
    if (!globalMasterIndexPtr_)
    {
        calcGlobalMasterIndex();
    }


    return *globalMasterIndexPtr_;
}

const Foam::List<Foam::word>&
Foam::solidGeneralContactFvPatchVectorField::shadowPatchNames() const
{
    if (!shadowPatchNamesPtr_)
    {
        calcShadowPatchNames();
    }

    return *shadowPatchNamesPtr_;
}

const Foam::List<Foam::label>&
Foam::solidGeneralContactFvPatchVectorField::shadowPatchIndices() const
{
	Info<<"In shadowPatchIndices() line "<<__LINE__<<endl;
    if (!shadowPatchIndicesPtr_)
    {
		Info<<"In shadowPatchIndices() line "<<__LINE__<<endl;
        calcShadowPatchNames();
    }
Info<<"*shadowPatchIndicesPtr_ in shadowPatchIndices() "<<*shadowPatchIndicesPtr_<<endl;
    return *shadowPatchIndicesPtr_;
}	

const Foam::List<Foam::word>&
Foam::solidGeneralContactFvPatchVectorField::shadowZoneNames() const
{
    if (!shadowZoneNamesPtr_)
    {
        calcShadowZoneNames();
    }

    return *shadowZoneNamesPtr_;
}

const Foam::List<Foam::label>&
Foam::solidGeneralContactFvPatchVectorField::shadowZoneIndices() const
{
    if (!shadowZoneIndicesPtr_)
    {
        calcShadowZoneNames();
    }

    return *shadowZoneIndicesPtr_;
}

// Map from self
void Foam::solidGeneralContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    FatalErrorIn
    (
        "void Foam::solidGeneralContactFvPatchVectorField::autoMap"
        "("
        "    const fvPatchFieldMapper& m"
        ")"
    )   << "member mapping not implemented" << endl;

    solidTractionFvPatchVectorField::autoMap(m);
}

// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::solidGeneralContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    FatalErrorIn
    (
        "void Foam::solidGeneralContactFvPatchVectorField::rmap"
        "("
        "    const fvPatchField<vector>& ptf,"
        "    const labelList& addr"
        ")"
    )   << "member mapping not implemented" << endl;

    solidTractionFvPatchVectorField::rmap(ptf, addr);
	// not sure if pointers are mapped correctly
    // be careful when there are topological changes to the patch
}

void solidGeneralContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

//**************************************************** START General**********************************************	
	Info<<"In updateCoeffs() line :"<<__LINE__<<endl;
//	Info<<"shadowPatchNames().size() in updateCoeffs():"<<shadowPatchNames().size()<<endl;
	boolList activeContactPairs(shadowPatchNames().size(), false);
	
	    // if it is a new time step then reset iCorr
    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Delete friction heat rate to force its recalculation when thermal
        // boundaries ask for it
        deleteDemandDrivenData(QcPtr_);
        deleteDemandDrivenData(QcsPtr_);

        if (globalMaster())
        {
			Info<<"Does it enter this check? in updateCoeffs():"<<__LINE__<<endl;
            forAll(activeContactPairs, shadowI)
            {
				Info<<"Which shadowI? "<<shadowI<<endl;
                // Let the contact models know that it is a new time-step, in
                // case they need to update anything
                normalModel(shadowI).newTimeStep();
                frictionModel(shadowI).newTimeStep();
				
		// **************** based on solid4Foam ****************
		// zoneToZonesNewGgi()[shadowI].clearPrevCandidateMasterNeighbors();
		// **************** end solid4Foam ****************
			}
        }
    }
	
	// Method
    // Move all global face zones to the deformed configuration
    // Clear interpolator weights
    // Perform quick check to find potential contacting pairs
    // Call normal and frction contact models for active contacting pairs
    // Accumulate contact force contributions for all active contact pairs


    if (rigidMaster_)
    {
        // Set to master to traction free to mimic a rigid patch
        traction() = vector::zero;
    }	
	else
    {
	// Move all global face zones to the deformed configuration
        if (globalMaster())
        {
			Info<<"Step1: Calling moveFaceZonesToDeformedConfiguration() in updateCoeffs():"<<__LINE__<<endl;
           // Move the master and slave zone to the deformed configuration
            moveFaceZonesToDeformedConfiguration();
			Info<<"Which patch was moved? "<<patch().name()<<endl;
        }


		// Clear interpolator weights
		Info<<"What is the size of activeContactPairs? "<<shadowPatchNames().size()<<endl;
        forAll(activeContactPairs, slaveI)
        {
			Info<<"Which slaveI? "<<slaveI<<endl;
            if (localSlave()[slaveI])
            {
				Info<<"In updateCoeffs():"<<__LINE__<<endl;
                zoneToZoneNewGgi(slaveI).movePoints
                (
                    tensorField(0), tensorField(0), vectorField(0)
                );
			Info<<"In updateCoeffs():"<<__LINE__<<endl;
            }
        }
		

		// Accumulated traction for the current patch
        vectorField curPatchTraction(patch().size(), vector::zero);
	//	Info<<"What is curPatchTraction? "<<curPatchTraction<<endl;
		Info<<"What is the size of the patch? "<<patch().size()<<endl;
		
		// Only the local masters calculates the contact force and the local
        // master interpolates this force
        const boolList& locSlave = localSlave();
		Info<<"Step1.1: Which patch is the local master? "<<patch().name()<<endl;
		Info<<"What is the locSlave boolList? "<<locSlave<<endl;
		
        // Create master bounding box used for quick check
        boundBox masterBb(zone().localPoints(), false);
		
//		Info<<"What are zone().localPoints()? "<<zone().localPoints()<<endl;

        // The BB may have zero thickness in one of the directions e.g. for a
        // flat patch, so we will check for this and, if found, create an offset
        const scalar bbOff = bbOffset();
		Info<<"What is bbOff? in updateCoeffs():"<<bbOff<<endl;
		if (masterBb.minDim() < bbOff)
        {
            const vector bbDiag = masterBb.max() - masterBb.min();

            if (bbDiag.x() < bbOff)
            {
				Info<<"Does it enter this check? in updateCoeffs():"<<__LINE__<<endl;
                vector offset(bbOff, 0, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
				Info<<"What is offset? in updateCoeffs():"<<offset<<endl;
            }
            else if (bbDiag.y() < bbOff)
            {
                vector offset(0, bbOff, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
            else if (bbDiag.z() < bbOff)
            {
                vector offset(0, 0, bbOff);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
        }
		
		forAll(activeContactPairs, shadowI)
        {
			Info<<"In updateCoeffs():"<<__LINE__<<endl;
			Info<<"What is activeContactPairs? in updateCoeffs():"<<activeContactPairs<<endl;
            // Perform quick check to find potential contacting pairs
            // The quick check is based on the bounding box (BB) of the contact
            // pairs: if the BBs of pair intersect then we will designate the
            // pair as active.
			
			// Create shadow bounding box
            boundBox shadowBb(shadowZone(shadowI).localPoints(), false);
			
			// Check for a zero dimension in the shadowBb
            if (shadowBb.minDim() < bbOff)
            {
                const vector bbDiag = shadowBb.max() - shadowBb.min();

                if (bbDiag.x() < bbOff)
                {
                    vector offset(bbOff, 0, 0);
                    shadowBb.min() -= offset;
                    shadowBb.max() += offset;
                }
                else if (bbDiag.y() < bbOff)
                {
                    vector offset(0, bbOff, 0);
                    shadowBb.min() -= offset;
                    shadowBb.max() += offset;
                }
                else if (bbDiag.z() < bbOff)
                {
                    vector offset(0, 0, bbOff);
                    shadowBb.min() -= offset;
                    shadowBb.max() += offset;
                }
            }
			
			if (masterBb.overlaps(shadowBb))
            {
				Info<<"In updateCoeffs():"<<__LINE__<<endl;
                activeContactPairs[shadowI] = true;
            }
			
			Info<<"activeContactPairs in updateCoeffs():"<<activeContactPairs<<endl;
			// Call normal and frction contact models for active contacting
            // pairs
            // Accumulate contact force contributions for all active contact
            // pairs

            if (activeContactPairs[shadowI])
            {
				Info<<"START activeContactPairs[shadowI] check in updateCoeffs():"<<__LINE__<<endl;
				if (locSlave[shadowI])
                {
					Info<<"SLAVE of LOCAL pair in updateCoeffs():"<<__LINE__<<endl;
					Info<<"locSlave[shadowI] check in updateCoeffs():"<<__LINE__<<endl;
                    Info<< "The current patch in updateCoeffs() is "<< patch().name()<< endl;			
		Info<< "The size of current patch in updateCoeffs() is "<< patch().size()<< endl;
					Info<<"What is the locSlave[shadowI] boolList? "<<locSlave[shadowI]<<endl;
					// Correct normal and friction contact models for the
                    // current contact pair

                    // Calculate the slave patch face unit normals as they are
                    // units by both the normal and friction models
					const vectorField shadowPatchFaceNormals =
                        patchField
                        (
                            shadowPatchIndices()[shadowI],
                            shadowZoneIndices()[shadowI],
                            shadowZone(shadowI).faceNormals()
                        );

                    // Interpolate the master displacement increment to the
                    // slave patch as it is required by specific normal and
                    // friction contact models

                    vectorField patchDD(patch().size(), vector::zero);
					Info<<"What is the patchDD? "<<patchDD<<endl;
                    vectorField shadowPatchDD
                    (
                        patch().boundaryMesh()
                        [
                            shadowPatchIndices()[shadowI]
                        ].size(),
                        vector::zero
                    );
					
					if (movingMesh())
                    {
                        // Updated Lagrangian, we will directly lookup the
                        // displacement increment

                        const volVectorField& DD =
                            db().lookupObject<volVectorField>("DU");

                        patchDD = DD.boundaryField()[patch().index()];
                        shadowPatchDD =
                            DD.boundaryField()[shadowPatchIndices()[shadowI]];
                    }
					else
                    {
                        // We will lookup the total displacement and old total
                        // displacement

                        const volVectorField& D =
                            db().lookupObject<volVectorField>("U");

                        patchDD =
                            D.boundaryField()[patch().index()]
                          - D.oldTime().boundaryField()[patch().index()];
                        shadowPatchDD =
                            D.boundaryField()[shadowPatchIndices()[shadowI]]
                          - D.oldTime().boundaryField()
                            [
                                shadowPatchIndices()[shadowI]
                            ];
                    }
					
		Info<<"SLAVE of LOCAL pair in updateCoeffs()"<<__LINE__<<endl;
			// Master zone DD
                    const vectorField zoneDD =
                        zoneField
                        (
                            zoneIndex(),
                            patch().index(),
						    patchDD
                        );
					
		Info<< "The current field in updateCoeffs() is "<< dimensionedInternalField().name()<< endl;
		Info<< "The current dimensionedInternalField().size() in updateCoeffs() is "<< dimensionedInternalField().size()<< endl;
		Info<< "The current patch in updateCoeffs() is "<< patch().name()<< endl;			
		Info<< "The size of current patch in updateCoeffs() is "<< patch().size()<< endl;
		
                    
		Info<<"updateCoeffs() - zoneDD.size(): "<<zoneDD.size()<<endl;
		Info<<"updateCoeffs() - zone().size(): "<<zone().size()<<endl;
	//	Info<<"updateCoeffs() - shadowZone(shadowI).size(): "<<shadowZone(shadowI).size()<<endl;
		Info<<"updateCoeffs() - patchDD.size(): "<<patchDD.size()<<endl;		
	//	Info<<"shadowPatchIndices()[shadowI]: "<<shadowPatchIndices()[shadowI]<<endl;
	//	Info<<"shadowZoneIndices()[shadowI]: "<<shadowZoneIndices()[shadowI]<<endl;
		Info<< "shadowI in updateCoeffs() "<<shadowI<< endl; 
		
		Info<<"SLAVE of LOCAL pair in updateCoeffs() line "<<__LINE__<<endl;	
                    // Master patch DD interpolated to the slave patch
                    
					const vectorField patchDDInterpToShadowPatch =
                        patchField
                        (
                            shadowPatchIndices()[shadowI],
                            shadowZoneIndices()[shadowI],
							//checking shadowZone(shadowI) instead of zoneToZone(shadowI)
                        //    zoneToZoneNewGgi(shadowI).masterToSlave(zoneDD)()
							zoneToZoneNewGgi(shadowI).slaveToMaster(zoneDD)()
                        );
										
					
					// *************** start ERROR (noMatchingFunctionCall)******************
					/*
					FatalError
                        << "Disabled: use jasakSolidContact" << abort(FatalError);
                     */ 
					Info<<"SLAVE of LOCAL pair in updateCoeffs() line "<<__LINE__<<endl;
					// normalModel()[shadowI].correct 
					 normalModel(shadowI).correct
                     (
                         shadowPatchFaceNormals,
                    //     zoneToZoneNewGgi(shadowI),
					/*	shadowZonesNewGgi()[shadowI].globalPointToPatch
									(
										zoneToZonesNewGgi()[shadowI].slavePointDistanceToIntersection()
									), */ 
					//	zoneToZoneNewGgi(shadowI).masterPointDistanceToIntersection(),
						zoneToZoneNewGgi(shadowI).slavePointDistanceToIntersection(),
                         shadowPatchDD,
                         patchDDInterpToShadowPatch
                     );
					 
					
					// *************** end ERROR (noMatchingFunctionCall)******************
					Info<<"SLAVE of LOCAL pair in updateCoeffs()"<<__LINE__<<endl;
					
					frictionModel(shadowI).correct
                    (
                        normalModel(shadowI).slavePressure(),
                        shadowPatchFaceNormals,
                        normalModel(shadowI).areaInContact(),
                        shadowPatchDD,
                        patchDDInterpToShadowPatch
                    );
					
					Info<<"SLAVE of LOCAL pair in updateCoeffs()"<<__LINE__<<endl;
					// Accumulate traction

                    curPatchTractions(shadowI) =
                        frictionModel(shadowI).slaveTraction()
                        + normalModel(shadowI).slavePressure();

                    curPatchTraction += curPatchTractions(shadowI);
				Info<<"curPatchTraction local Slave in updateCoeffs()"<<curPatchTraction<<endl;
				}

				else // local master
                {
					Info<<"LOCAL pair MASTER in updateCoeffs()"<<__LINE__<<endl;
					Info<< "The current patch in updateCoeffs() is "<< patch().name()<< endl;			
					Info<< "The size of current patch in updateCoeffs() is "<< patch().size()<< endl;
					Info<< "The current field in updateCoeffs() is "<<dimensionedInternalField().name()<< endl;
					// Get traction from local slave

                    const volVectorField& field =
                        db().lookupObject<volVectorField>
                        (
                            dimensionedInternalField().name()
                        );
					
					Info<<"LOCAL pair MASTER in updateCoeffs()"<<__LINE__<<endl;
						
					const solidGeneralContactFvPatchVectorField&
                        localMasterField =
                        refCast<const solidGeneralContactFvPatchVectorField>
                        (
                            field.boundaryField()
                            [
                                shadowPatchIndices()[shadowI]
                            ]
                        );
						
					Info<<"LOCAL pair MASTER in updateCoeffs()"<<__LINE__<<endl;
				//	Info<<"What is field.boundaryField()? in updateCoeffs()"<<field.boundaryField()[shadowPatchIndices()[shadowI]]<<endl;
					Info<<"*localMasterField in updateCoeffs()"<<*localMasterField<<endl;
					Info<<"shadowI in updateCoeffs()"<<shadowI<<endl;
				//	Info<<"findShadowID(patch().index()) in updateCoeffs()"<<findShadowID(patch().index())<<endl;
					
					const label masterShadowI =
                        localMasterField.findShadowID(patch().index());
					
					Info<<"masterShadowI in updateCoeffs()"<<masterShadowI<<endl;
                
					vectorField shadowPatchTraction =
                        -localMasterField.frictionModel
                        (
                            masterShadowI
                        ).slaveTractionForMaster()   //slaveTractionForMaster()
                        -localMasterField.normalModel
                        (
                            masterShadowI
                        ).slavePressure();
						
					Info<<"LOCAL pair MASTER in updateCoeffs()"<<__LINE__<<endl;
					
					//************** START (remove this evaluation)*******
					
					vectorField shadowZoneTraction =
                        zoneField
                        (
                            shadowZoneIndices()[shadowI],
                            shadowPatchIndices()[shadowI],
                            shadowPatchTraction
                        ); 
						
					
					Info<<"LOCAL pair MASTER in updateCoeffs()"<<__LINE__<<endl;
					Info<<"size of shadowZoneTraction? in updateCoeffs()"<<shadowZoneTraction.size()<<endl;
										
					// Face-to-face
					vectorField masterZoneTraction =
                        localMasterField.zoneToZoneNewGgi
                        (
                            masterShadowI
                        ).masterToSlave(shadowPatchTraction)();
					//	).slaveToMaster(shadowZoneTraction)();	//.slavePointDistanceToIntersection()
						
					Info<<"LOCAL pair MASTER in updateCoeffs()"<<__LINE__<<endl;
					Info<<"shadowI in updateCoeffs()"<<shadowI<<endl;
					Info<<"masterShadowI in updateCoeffs()"<<masterShadowI<<endl;
					Info<<"What is masterZoneTraction? in updateCoeffs()"<<masterZoneTraction<<endl;
					Info<<"size of masterZoneTraction? in updateCoeffs()"<<masterZoneTraction.size()<<endl;
					Info<<"LOCAL pair MASTER in updateCoeffs()"<<__LINE__<<endl;
					
                    // We store master patch traction as thermalGeneralContact
                    // uses it
					
                    curPatchTractions(shadowI) = 	
                        patchField
                        (
                            patch().index(),
                            zoneIndex(),
                            masterZoneTraction
                        ); 
							
					// ********* END (remove this evaluation)***********
					
					Info<<"LOCAL pair MASTER in updateCoeffs()"<<__LINE__<<endl;	
				//	Info<<"curPatchTractions(shadowI) in updateCoeffs()"<<curPatchTractions(shadowI)<<endl;
                    curPatchTraction += curPatchTractions(shadowI);
					Info<<"curPatchTraction local Master in updateCoeffs()"<<curPatchTraction<<endl;
				}				
			} // if contact pair is active
		} // forAll contact pairs
		
		// Set master gradient based on accumulated traction
        traction() = curPatchTraction;
	}
	
//**************************************************** END General**********************************************
       Info<<"In updateCoeffs() line "<<__LINE__<<endl;
	   	Info<<"traction() in updateCoeffs() 2045 "<<traction()<<endl;

    solidTractionFvPatchVectorField::updateCoeffs();
	Info<<"In updateCoeffs() line "<<__LINE__<<endl;
}

const Foam::scalarField& Foam::solidGeneralContactFvPatchVectorField::Qc() const
{
    if (!QcPtr_)
    {
        calcQc();
    }

    return *QcPtr_;
}

//- Increment of dissipated energy due to friction for each pair
const Foam::scalarField& Foam::solidGeneralContactFvPatchVectorField::Qc
(
    const label shadowI
) const
{
    if (!QcsPtr_)
    {
        calcQcs();
    }

    return (*QcsPtr_)[shadowI];
}

void Foam::solidGeneralContactFvPatchVectorField::calcQcs() const
{
    const boolList& locSlave = localSlave();

    QcsPtr_ = new List<scalarField>(locSlave.size());

    bool sigmaCauchyFound =
        db().foundObject<volSymmTensorField>
        (
            "sigmaCauchy"
        );

    const scalar deltaT = patch().boundaryMesh().mesh().time().deltaT().value();

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    forAll(locSlave, shadowI)
    {
        scalarField& Qc = (*QcsPtr_)[shadowI];
        Qc.setSize(patch().size(), 0.0);

        // For now, we assume traction is constant over time-step
        // Todo: use trapezoidal rule
        vectorField curTraction(Qc.size(), vector::zero);

        // sigma/sigmaCauchy is up-to-date as Qc is called after momentum loop
        // has converged and sigma has been updated and mesh moved
        if (sigmaCauchyFound)
        {
            const symmTensorField& sigma =
                db().lookupObject<volSymmTensorField>
                (
                    "sigmaCauchy"
                ).boundaryField()[patch().index()];

            curTraction = patch().nf() & sigma;
        }
        else
        {
            const symmTensorField& sigma =
            db().lookupObject<volSymmTensorField>
            (
                "sigma"
            ).boundaryField()[patch().index()];

            curTraction = patch().nf() & sigma;
        }

        // Calculate Qc for shadowI

        vectorField curPatchSlip(Qc.size(), vector::zero);

        // Calculate slip
        if (locSlave[shadowI])
        {
			
            curPatchSlip = frictionModel(shadowI).slip();
			
        }
        else
        {
            const solidGeneralContactFvPatchVectorField& shadowPatchField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadowI]]
            );

            const label locShadowID =
            shadowPatchField.findShadowID(patch().index());

                        
			vectorField shadowPatchSlip =
                shadowPatchField.frictionModel(locShadowID).slip();
				
			
            
			vectorField shadowZoneSlip =
                zoneField
                (
                    shadowZoneIndices()[shadowI],
                    shadowPatchIndices()[shadowI],
                    shadowPatchSlip
                );
			

            // Interpolate from shadow to the current patch
            // Face-to-face
			
            
			vectorField curZoneSlip =
                shadowPatchField.zoneToZoneNewGgi(locShadowID).slaveToMaster
                (
                    shadowZoneSlip
                );
			
						
				
			
            
			curPatchSlip =
            patchField
            (
                patch().index(),
                zoneIndex(),
                curZoneSlip
            );
				
			
        }

        // Heat flux rate: rate of dissipated frictional energy
        // The dot product of the traction vectors and the slip vectors
        // gives the dissipated frictional energy rate per unit area, which
        // is always positive
        Qc = mag(curTraction & (curPatchSlip/deltaT));
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcBbOffset() const
{
    if (bbOffset_ != 0)
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::solidGeneralContactFvPatchVectorField::"
            "calcBbOffset() const"
        )   << "already set" << abort(FatalError);
    }

    // We will set the BB offset to five times the average dimension of the
    // smallest face on the zone

    scalar minDim = GREAT;

    if (patch().size() > 0)
    {
        minDim = min(sqrt(patch().magSf()));
    }

    bbOffset_ = 5.0*returnReduce(minDim, minOp<scalar>());

    if (debug)
    {
        Info<< nl << "The bbOffset is " << bbOffset_ << endl;
    }
}

Foam::scalar Foam::solidGeneralContactFvPatchVectorField::bbOffset() const
{
    if (bbOffset_ == 0)
    {
        calcBbOffset();
    }

    return bbOffset_;
}


const Foam::vectorField&
Foam::solidGeneralContactFvPatchVectorField::curPatchTractions
(
    const label shadowI
) const
{
	Info<<"In curPatchTractions(..) line: "<<__LINE__<<endl;
    if (!curPatchTractionPtr_)
    {
        makeCurPatchTractions();
    }

    return (*curPatchTractionPtr_)[shadowI];
}


Foam::vectorField&
Foam::solidGeneralContactFvPatchVectorField::curPatchTractions
(
    const label shadowI
)
{
	Info<<"In curPatchTractions(..) line: "<<__LINE__<<endl;
    if (!curPatchTractionPtr_)
    {
		Info<<"In curPatchTractions(..) line: "<<__LINE__<<endl;
        makeCurPatchTractions();
    }
	
	Info<<"In curPatchTractions(..) line: "<<__LINE__<<endl;
	Info<< "The current patch in curPatchTractions() is "<< patch().name()<< endl;			
	Info<< "patch().size() in curPatchTractions() is "<< patch().size()<< endl;
	
	Info<<"(*curPatchTractionPtr_) in curPatchTractions(..) "<<(*curPatchTractionPtr_)<<endl;
    return (*curPatchTractionPtr_)[shadowI];
}



void Foam::solidGeneralContactFvPatchVectorField::calcZoneToZones() const
{
	Info<< "The current field in calcZoneToZones() is "<< dimensionedInternalField().name()<< endl;
	Info<< "The current patch in calcZoneToZones() is "<< patch().name()<< endl;			
	Info<< "The size of current patch in calcZoneToZones() is "<< patch().size()<< endl;
    // Create zone-to-zone interpolation
    if (!zoneToZonesNewGgi_.empty())
    {
        FatalErrorIn
        (
            "void solidGeneralContactFvPatchVectorField::calcZoneToZones()"
            "const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    zoneToZonesNewGgi_.setSize(shadowPatchNames().size());

    const boolList& locSlave = localSlave();
	
	Info<<"What does the locSlave boolList return? "<<locSlave<<endl;

    forAll(zoneToZonesNewGgi_, shadowI)
    {
        // Only the local slave creates the interpolator
        if (locSlave[shadowI])
        {
			Info<<"locSlave[shadowI] check in calcZoneToZones() line: "<<__LINE__<<endl;
            Info<<"What is the locSlave[shadowI] boolList? "<<locSlave[shadowI]<<endl;
			
            zoneToZonesNewGgi_.set
                (
                    shadowI,
                    new newGgiStandAlonePatchInterpolation 
                    (
                        shadowZone(shadowI), // master
                        zone(), // slave
					//	shadowZone(shadowI), // master
                        tensorField(0),
                        tensorField(0),
                        vectorField(0), // Slave-to-master separation.
                        true,           // global data
                        0,              // Non-overlapping face tolerances
                        0,              //
                        true,           // Rescale weighting factors.
                        newGgiInterpolation::AABB   
                        //newGgiInterpolation::BB_OCTREE
                        //newGgiInterpolation::THREE_D_DISTANCE
                        //newGgiInterpolation::N_SQUARED
                    )
                );
			Info<<"In calcZoneToZones() line: "<<__LINE__<<endl;	
        }
    }
}


// ******************* Definition from solid4foam***********************************************

const Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidGeneralContactFvPatchVectorField::shadowZonesNewGgi() const
{
   Info<<"CHECK solid4Foam  Here I am in shadowZonesNewGgi()"<<__LINE__<<endl;
   if (globalMasterPtr_)
    {
		Info<<"Here I am in shadowZonesNewGgi()"<<__LINE__<<endl;
        if (shadowZonesNewGgi_.empty())
        {
            calcShadowZonesNewGgi();
        }

        return shadowZonesNewGgi_;
    }
    else
    {
		Info<<"Here I am in shadowZonesNewGgi()"<<__LINE__<<endl;
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        const solidGeneralContactFvPatchVectorField& shadowPatchField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );

        return shadowPatchField.shadowZonesNewGgi();
    }
	
}


Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidGeneralContactFvPatchVectorField::shadowZonesNewGgi()
{
	Info<<"CHECK solid4Foam Here I am in shadowZonesNewGgi()"<<__LINE__<<endl;
    if (globalMasterPtr_)
    {
		Info<<"Here I am in shadowZonesNewGgi()"<<__LINE__<<endl;
        if (shadowZonesNewGgi_.empty())
        {
            calcShadowZonesNewGgi();
        }

        return shadowZonesNewGgi_;
    }
    else
    {
		Info<<"Here I am in shadowZonesNewGgi()"<<__LINE__<<endl;
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        solidGeneralContactFvPatchVectorField& shadowPatchField =
            const_cast<solidGeneralContactFvPatchVectorField&>
            (
                refCast<const solidGeneralContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.shadowZonesNewGgi();
    }
}

// ******************* End definition from solid4foam***********************************************

// ******************* Definition from solid4foam***********************************************
/* const Foam::globalPolyPatch&
Foam::solid4ContactFvPatchVectorField::zoneNewGgi() const
{

}


Foam::globalPolyPatch& Foam::solid4ContactFvPatchVectorField::zoneNewGgi()
{

}
*/
// ******************* End definition from solid4foam***********************************************


// ******************* Definition from solid4foam***********************************************
const Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
Foam::solidGeneralContactFvPatchVectorField::zoneToZonesNewGgi() const
{
	Info<<"CHECK solid4Foam Here I am in solid4Foam's zoneToZonesNewGgi()"<<__LINE__<<endl;
    if (master_)
    {
		Info<<"CHECK solid4Foam Here I am in solid4Foam's zoneToZonesNewGgi()"<<__LINE__<<endl;
        if (zoneToZonesNewGgi_.empty())
        {
			Info<<"CHECK solid4Foam Here I am in solid4Foam's zoneToZonesNewGgi()"<<__LINE__<<endl;
            calcZoneToZones();
        }
		
	Info<< "CHECK The current patch is "<< patch().name()<< endl;	
        return zoneToZonesNewGgi_;    // This needs to be FIXED later 
    }
    else
    {
		Info<<"CHECK solid4Foam Here I am in solid4Foam's zoneToZonesNewGgi()"<<__LINE__<<endl;
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        const solidGeneralContactFvPatchVectorField& shadowPatchField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );
			
	Info<< "CHECK The current patch is "<< patch().name()<< endl;	
        return shadowPatchField.zoneToZonesNewGgi();
    }
}


Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
Foam::solidGeneralContactFvPatchVectorField::zoneToZonesNewGgi()
{
	Info<<"CHECK solid4Foam Here I am in solid4Foam's zoneToZonesNewGgi()"<<__LINE__<<endl;
    if (globalMasterPtr_)
    {
		Info<<"CHECK solid4Foam Here I am in solid4Foam's zoneToZonesNewGgi()"<<__LINE__<<endl;
        if (zoneToZonesNewGgi_.empty())
        {
			Info<<"Does it enter solid4Foam's zoneToZonesNewGgi()? CHECK line "<<__LINE__<<endl;
            calcZoneToZones();
        }
		Info<<"zoneToZonesNewGgi_ in solid4Foam's zoneToZonesNewGgi() line "<<__LINE__<<endl;
	
         return zoneToZonesNewGgi_;    // This needs to be FIXED later 
    }
    else
    {
		Info<<"CHECK solid4Foam Here I am in zoneToZonesNewGgi()"<<__LINE__<<endl;
        // We will const_cast the shadow patch so we can delete the weights when
        // the zones move
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        solidGeneralContactFvPatchVectorField& shadowPatchField =
            const_cast<solidGeneralContactFvPatchVectorField&>
            (
                refCast<const solidGeneralContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.zoneToZonesNewGgi();
    }

}

// ****************************** end definition from solid4foam*********************************

const Foam::newGgiStandAlonePatchInterpolation&
Foam::solidGeneralContactFvPatchVectorField::zoneToZoneNewGgi
(
    const label shadowI
) const
{ 
	Info<<"In zoneToZoneNewGgi()"<<__LINE__<<endl;
   if (!localSlave()[shadowI])
    {
        FatalErrorIn("zoneToZoneNewGgi(const label shadowI)")
            << "Only the local slave can call the zoneToZoneNewGgi interpolator"
            << abort(FatalError);
    }
	Info<<"In zoneToZoneNewGgi()"<<__LINE__<<endl;
    if (zoneToZonesNewGgi_.empty())
    {
        word zoneName =
            patch().boundaryMesh().mesh().faceZones()[zoneIndex()].name();

        if (debug)
        {
            Info<< "Initializing the GGI interpolators for " << zoneName
                << endl;
        }

        calcZoneToZones();
    }
	Info<<"In zoneToZoneNewGgi()"<<__LINE__<<endl;

    return zoneToZonesNewGgi_[shadowI];
}

Foam::newGgiStandAlonePatchInterpolation&
Foam::solidGeneralContactFvPatchVectorField::zoneToZoneNewGgi(const label shadowI)
{
	Info<<"In zoneToZoneNewGgi()"<<__LINE__<<endl;
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("zoneToZone(const label shadowI)")
            << "Only the local slave can call the zoneToZone interpolator"
            << abort(FatalError);
    }
	
	Info<<"In zoneToZoneNewGgi()"<<__LINE__<<endl;

    if (zoneToZonesNewGgi_.empty())
    {
        word zoneName =
            patch().boundaryMesh().mesh().faceZones()[zoneIndex()].name();

        if (debug)
        {
            Info<< "Initializing the GGI interpolators for " << zoneName
                << endl;
        }

        calcZoneToZones();
    }
	
	Info<<"In zoneToZoneNewGgi()"<<__LINE__<<endl;
//	Info<<"zoneToZonesNewGgi_[shadowI] in zoneToZoneNewGgi()"<<zoneToZonesNewGgi_<<endl;

    return zoneToZonesNewGgi_[shadowI];
}



bool Foam::solidGeneralContactFvPatchVectorField::globalMaster() const
{
    if (!globalMasterPtr_)
    {
        calcGlobalMaster();
    }

	Info<< "The current field in globalMaster() is "<< dimensionedInternalField().name()<< endl;
	Info<< "The current patch in globalMaster()) is "<< patch().name()<< endl;			
	Info<< "The size of current patch in globalMaster() is "<< patch().size()<< endl;
	Info<< "*globalMasterPtr_ in globalMaster() is "<<*globalMasterPtr_<< endl;
    return *globalMasterPtr_;
}


/*
//  Move the contact face zone patches to the deformed position
void solidGeneralContactFvPatchVectorField::moveFaceZonePatches()
{
   if (slaveToMasterPatchToPatchInterpolatorPtr_)
  {
      delete slaveToMasterPatchToPatchInterpolatorPtr_;
      slaveToMasterPatchToPatchInterpolatorPtr_ =
          new PatchToPatchInterpolation<
              PrimitivePatch<
                  face, Foam::List, pointField
                  >, PrimitivePatch<face, Foam::List, pointField> >
          (
              *slaveFaceZonePatchPtr_, // from zone
              *masterFaceZonePatchPtr_, // to zone
              alg_,
              dir_
          );
  }
  else if (slaveToMasterGgiInterpolatorPtr_)
  {
      delete slaveToMasterGgiInterpolatorPtr_;
      slaveToMasterGgiInterpolatorPtr_ =
        new GGIInterpolation<
            PrimitivePatch<
                face, Foam::List, pointField
                >, PrimitivePatch< face, Foam::List, pointField > >
          (
              *masterFaceZonePatchPtr_, // master zone
              *slaveFaceZonePatchPtr_, // slave zone
              tensorField(0),
              tensorField(0),
              vectorField(0),
              0.0,
              0.0,
              true,
              ggiInterpolation::AABB
          );
  }

}
*/




// Write
void solidGeneralContactFvPatchVectorField::write(Ostream& os) const
{
  Info << "solidGeneralContactFvPatchVectorField::write..." << endl;
  
  //****************************START general***************************************//
  
  solidTractionFvPatchVectorField::write(os);  

	if(
	    !localSlavePtr_
		&& dimensionedInternalField().name() == "U_0"
		)
		{
//	  	Info<<"Here I am in first U_0 check in write()"<<__LINE__<<endl;
		return;
		} 

   os.writeKeyword("rigidMaster")
        << rigidMaster_ << token::END_STATEMENT << nl;

    // Write the dict from the first contact model

    const label shadowI = 0;

    if(!localSlavePtr_) //remove this check later, since localSlave should re-compute the local slave
        FatalError  << "solidGeneralContactFvPatchVectorField::write: localSlavePtr_ NOT defined:" 
                    << "Cannot write slave information because no slave identified!"  
                    << exit(FatalError);; 

    if (localSlave()[shadowI])
    {
        os.writeKeyword("generalNormalContactModel")
            << normalModel(shadowI).type() << token::END_STATEMENT << nl;
        normalModel(shadowI).writeDict(os);

        os.writeKeyword("generalFrictionContactModel")
            << frictionModel(shadowI).type() << token::END_STATEMENT << nl;
        frictionModel(shadowI).writeDict(os);
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                dimensionedInternalField().name()
            );

        const solidGeneralContactFvPatchVectorField& localSlaveField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()
                [
                    shadowPatchIndices()[shadowI]
                ]
            );

        const label localSlaveID =
            localSlaveField.findShadowID(patch().index());

        os.writeKeyword("generalNormalContactModel")
            << localSlaveField.normalModel(localSlaveID).type()
            << token::END_STATEMENT << nl;
        localSlaveField.normalModel(localSlaveID).writeDict(os);

        os.writeKeyword("generalFrictionContactModel")
            << localSlaveField.frictionModel(localSlaveID).type()
            << token::END_STATEMENT << nl;
        localSlaveField.frictionModel(localSlaveID).writeDict(os);
    }
	//****************************END general***************************************//
	
    
}

Foam::label Foam::solidGeneralContactFvPatchVectorField::findShadowID
(
    const label patchID
) const
{
    label shadowI = -1;

    const labelList shadowIDs = shadowPatchIndices();

    forAll(shadowIDs, I)
    {
        if (patchID == shadowIDs[I])
        {
            shadowI = I;
            break;
        }
    }

    if (shadowI == -1)
    {
        FatalErrorIn("findShadowID(const label patchID)")
            << "shadow patch not found!" << abort(FatalError);
    }

    return shadowI;
}


void Foam::solidGeneralContactFvPatchVectorField::makeCurPatchTractions() const
{
	Info<<"In makeCurPatchTractions() line "<<__LINE__<<endl;
	
	if (curPatchTractionPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "makeCurPatchTractions() const"
        )   << "curPatchTractionPtr_ already set" << abort(FatalError);
    }
	
	Info<< "The current patch in makeCurPatchTractions() is "<< patch().name()<< endl;			
	Info<< "vectorField(patch().size() in makeCurPatchTractions() is "<< patch().size()<< endl;
    Info<< "shadowPatchNames().size() in makeCurPatchTractions() is "<< shadowPatchNames().size()<< endl;

    curPatchTractionPtr_ =
        new List<vectorField>
        (
            shadowPatchNames().size(),
            vectorField(patch().size(), vector::zero)
        );
	Info<< "curPatchTractionPtr_ in makeCurPatchTractions() is "<<*curPatchTractionPtr_<< endl;
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
makePatchTypeField(fvPatchVectorField, solidGeneralContactFvPatchVectorField)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
