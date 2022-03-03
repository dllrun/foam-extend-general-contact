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
    solid4GeneralContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solid4GeneralContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"

#define ISDEBUG false
#define ActivePairDEBUG false
#define currentMasterDEBUG false
#define localSlaveDEBUG false
#define normalModelDEBUG false
#define masterOfPairDEBUG false
#define slaveOfPairDEBUG false
#define bbOffDEBUG false



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::solid4GeneralContactFvPatchVectorField::movingMesh() const
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


void Foam::solid4GeneralContactFvPatchVectorField::
moveZonesToDeformedConfiguration()
{
    // Only the master moves the zones
    if (!firstPatchInList())
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

    // Target patch and zone indices
    const labelList& activePatchIndices = targetPatchIndices();
    const labelList& shadZoneIndices = targetGlobalPolyPatchFaceZoneIndices();

    forAll(activePatchIndices, targetI)
    {
        // Assemble the zone face displacement field to move the zones
        vectorField zoneD(globalPolyPatchFaceZone().globalPatch().size(), vector::zero);
        vectorField targetGlobalPolyPatchFaceZoneD(targetGlobalPolyPatchFaceZone(targetI).globalPatch().size(), vector::zero);

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
            const vectorField& targetPatchDD =
                DD.boundaryField()[activePatchIndices[targetI]];

            zoneD = globalPolyPatchFaceZone().patchFaceToGlobal(patchDD);		
            //    zoneField(globalPolyPatchFaceZoneIndex(), patch().index(), patchDD);
			
            targetGlobalPolyPatchFaceZoneD = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(targetPatchDD);
			/*
                    zoneField
                    (
                        shadZoneIndices[targetI],
                        activePatchIndices[targetI],
                        targetPatchDD
                    );
			*/
        }
        else
        {
            // Non-moving mesh: we will move the zones by the total displacement

            // Lookup the current total displacement field
            const volVectorField& D = db().lookupObject<volVectorField>("U");

            // Take a reference to the patch face total displacement field
            const vectorField& patchD =
                D.boundaryField()[patch().index()];

            const vectorField& targetPatchDD =
                D.boundaryField()[activePatchIndices[targetI]];

            zoneD = globalPolyPatchFaceZone().patchFaceToGlobal(patchD);
                //zoneField(globalPolyPatchFaceZoneIndex(), patch().index(), patchD);
            targetGlobalPolyPatchFaceZoneD = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(targetPatchDD);
                /*
                zoneField
                (
                    shadZoneIndices[targetI],
                    activePatchIndices[targetI],
                    targetPatchDD
                );
				*/
        }

        // Interpolate the zone face field to the zone points
        const pointField zonePointD = globalPolyPatchFaceZone().interpolator().faceToPointInterpolate(zoneD);
            //zoneFaceToPointInterpolate(globalPolyPatchFaceZoneIndex(), zoneD, -1);
        const pointField targetGlobalPolyPatchFaceZonePointD =
			targetGlobalPolyPatchFaceZone(targetI).interpolator().faceToPointInterpolate
            (
              targetGlobalPolyPatchFaceZoneD
            );
			/*
			zoneFaceToPointInterpolate
            (
                shadZoneIndices[targetI], targetGlobalPolyPatchFaceZoneD, targetI
            );
			*/

        // The zone deformed points are the initial position plus the
        // displacement
        const pointField zoneNewPoints =
            mesh.faceZones()[globalPolyPatchFaceZoneIndex()]().localPoints()
          + zonePointD;
        const pointField targetGlobalPolyPatchFaceZoneNewPoints =
            mesh.faceZones()[shadZoneIndices[targetI]]().localPoints()
          + targetGlobalPolyPatchFaceZonePointD;

        // Move the zones

        // Remove zones weights
        if (targetI == 0)
        {
            globalPolyPatchFaceZone().movePoints(zoneNewPoints);
        }
        targetGlobalPolyPatchFaceZone(targetI).movePoints(targetGlobalPolyPatchFaceZoneNewPoints);

        // We need to use const_cast to move the standAlonePatch points as the
        // movePoints function only clears weights
        // Also, be careful to move the points are opposed to the localPoints
        if (targetI == 0)
        {
            const_cast<pointField&>(globalPolyPatchFaceZone().globalPatch().points()) = zoneNewPoints;
        }
        const_cast<pointField&>(targetGlobalPolyPatchFaceZone(targetI).globalPatch().points()) =
            targetGlobalPolyPatchFaceZoneNewPoints;
    }
}


void Foam::solid4GeneralContactFvPatchVectorField::calcFirstPatchInList() const
{
    if (firstPatchInListPtr_)
    {
        FatalErrorIn
            (
                "void Foam::solid4GeneralContactFvPatchVectorField::"
                "calcFirstPatchInList() const"
            )   << "firstPatchInListPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solid4GeneralContact patch i.e. the one
    // with the lowest patch index

    if (firstPatchInListIndex() == patch().index())
    {
        firstPatchInListPtr_ = new bool(true);
    }
    else
    {
        firstPatchInListPtr_ = new bool(false);
    }
}


void Foam::solid4GeneralContactFvPatchVectorField::calcFirstPatchInListIndex() const
{
    if (firstPatchInListIndexPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solid4GeneralContactFvPatchVectorField::"
            "calcFirstPatchInListIndex() const"
        )   << "firstPatchInListIndexPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solid4GeneralContact patch i.e. the one
    // with the lowest patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    firstPatchInListIndexPtr_ = new label(-1);
    label& gMasterID = *firstPatchInListIndexPtr_;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solid4GeneralContactFvPatchVectorField::typeName
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
            "solid4GeneralContactFvPatchVectorField::"
            "calcFirstPatchInListIndex() const"
        )   << "There are multiple global masters" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< nl << "The global master contact patch is "
            << patch().boundaryMesh()[gMasterID].name() << endl;
    }
}


void Foam::solid4GeneralContactFvPatchVectorField::calcLocalSlave() const
{
    if (localSlavePtr_)
    {
        FatalErrorIn
        (
            "label Foam::solid4GeneralContactFvPatchVectorField::"
            "calcLocalSlave() const"
        )   << "localSlavePtr_ already set" << abort(FatalError);
    }

    localSlavePtr_ = new boolList(targetPatchNames().size(), false);

    boolList& localSlave = *localSlavePtr_;

    forAll(localSlave, targetI)
    {
        if (patch().index() < targetPatchIndices()[targetI])
        {
            localSlave[targetI] = true;

            Info<< "solid4GeneralContact: "
                << targetPatchNames()[targetI] << " (master)" << " to "
                << patch().name() << " (slave)" << endl;
        }
    }
}


const Foam::boolList&
Foam::solid4GeneralContactFvPatchVectorField::localSlave() const
{
    if (!localSlavePtr_)
    {
        calcLocalSlave();
    }

    return *localSlavePtr_;
}


void Foam::solid4GeneralContactFvPatchVectorField::calctargetPatchNames() const
{
    if (targetPatchNamesPtr_ || targetPatchIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solid4GeneralContactFvPatchVectorField::"
            "calctargetPatchNames() const"
        )   << "targetPatchNames_ or targetPatchIndices_ already set"
            << abort(FatalError);
    }

    // Add each solid4GeneralContact patch in the order of increasing patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    // Count target patches

    label nShadPatches = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solid4GeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
            nShadPatches++;
        }
    }

    targetPatchNamesPtr_ = new wordList(nShadPatches);
    wordList& targetPatchNames = *targetPatchNamesPtr_;

    targetPatchIndicesPtr_ = new labelList(nShadPatches);
    labelList& targetPatchIndices = *targetPatchIndicesPtr_;

    // Record target patch names

    label targetI = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solid4GeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
            targetPatchNames[targetI] = patch().boundaryMesh()[patchI].name();

            targetPatchIndices[targetI++] = patchI;
        }
    }
}


void Foam::solid4GeneralContactFvPatchVectorField::calcTargetGlobalPolyPatchFaceZoneNames() const
{
    if (targetGlobalPolyPatchFaceZoneNamesPtr_ || targetGlobalPolyPatchFaceZoneIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solid4GeneralContactFvPatchVectorField::"
            "calcTargetGlobalPolyPatchFaceZoneNames() const"
        )   << "targetGlobalPolyPatchFaceZoneNames_ or targetGlobalPolyPatchFaceZoneIndices_ already set"
            << abort(FatalError);
    }

    const wordList& shadNames = targetPatchNames();

    targetGlobalPolyPatchFaceZoneNamesPtr_ = new wordList(shadNames.size());
    wordList& targetGlobalPolyPatchFaceZoneNames = *targetGlobalPolyPatchFaceZoneNamesPtr_;

    targetGlobalPolyPatchFaceZoneIndicesPtr_ = new labelList(shadNames.size());
    labelList& targetGlobalPolyPatchFaceZoneIndices = *targetGlobalPolyPatchFaceZoneIndicesPtr_;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    forAll(shadNames, targetI)
    {
        word zoneName = shadNames[targetI] + "FaceZone";

        faceZoneID zone(zoneName, mesh.faceZones());

        if (!zone.active())
        {
            FatalErrorIn("solid4GeneralContactFvPatchVectorField")
                << "Face zone name " << zoneName
                << " not found.  Please check your zone definition."
                << abort(FatalError);
        }

        targetGlobalPolyPatchFaceZoneNames[targetI] = zoneName;

        targetGlobalPolyPatchFaceZoneIndices[targetI] = zone.index();
    }
}


void Foam::solid4GeneralContactFvPatchVectorField::calcNormalModels(const dictionary& dict) const
{
    if (!normalModels_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::"
            "calcNormalModel() const"
        )   << "normalModels already set" << abort(FatalError);
    }

    normalModels_.setSize(targetPatchNames().size());

    const boolList& locSlave = localSlave();

    forAll(normalModels_, targetI)
    {
        // Only the local slave creates the contact model
        if (locSlave[targetI])
        {
			const dictionary* contactDictPtr = NULL;
			if (dict.found("generalNormalContactModel") )
            {				
                contactDictPtr = &dict;
				
			}
			
			const dictionary& contactDict = *contactDictPtr;
			
            // Calculate normal contact forces
            normalModels_.set
            (
                targetI,
                generalNormalContactModel::New
                (
                    word(contactDict.lookup("generalNormalContactModel")),
                    patch().boundaryMesh()[targetPatchIndices()[targetI]],
                    contactDict,
                    targetPatchIndices()[targetI], // master
                    patch().index(), // slave
                    targetGlobalPolyPatchFaceZone(targetI).globalPatch(), // master
                    globalPolyPatchFaceZone().globalPatch() // slave
                )
            );
        }
    }
}


Foam::generalNormalContactModel&
Foam::solid4GeneralContactFvPatchVectorField::normalModel(const label targetI)
{
    if (!localSlave()[targetI])
    {
        FatalErrorIn("normalModel(const label targetI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (normalModels_.empty())
    {
        calcNormalModels(dict_);
    }

    return normalModels_[targetI];
}


const Foam::generalNormalContactModel&
Foam::solid4GeneralContactFvPatchVectorField::normalModel
(
    const label targetI
) const
{
    if (!localSlave()[targetI])
    {
        FatalErrorIn("normalModel(const label targetI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (normalModels_.empty())
    {
        calcNormalModels(dict_);
    }

    return normalModels_[targetI];
}


void Foam::solid4GeneralContactFvPatchVectorField::calcFrictionModels(const dictionary& dict) const
{
    if (!frictionModels_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::"
            "calcFrictionModel(targetI) const"
        )   << "frictionModelPtr_[targetI] already set" << abort(FatalError);
    }

    frictionModels_.setSize(targetPatchNames().size());

    const boolList& locSlave = localSlave();

    forAll(frictionModels_, targetI)
    {
        if (locSlave[targetI])
        {
			const dictionary* contactDictPtr = NULL;
			if (dict.found("generalFrictionContactModel") )
            {				
                contactDictPtr = &dict;
				
			}
			
			const dictionary& contactDict = *contactDictPtr;
			
            frictionModels_.set
                (
                    targetI,
                    generalFrictionContactModel::New
                    (
                        word(contactDict.lookup("generalFrictionContactModel")),
                        patch().boundaryMesh()[targetPatchIndices()[targetI]],
                        contactDict,
                        targetPatchIndices()[targetI], // master
                        patch().index() // slave
                    )
                );
        }
    }
}


Foam::generalFrictionContactModel&
Foam::solid4GeneralContactFvPatchVectorField::frictionModel(const label targetI)
{
    if (!localSlave()[targetI])
    {
        FatalErrorIn("frictionModel(const label targetI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (frictionModels_.empty())
    {
        calcFrictionModels(dict_);
    }

    return frictionModels_[targetI];
}


const Foam::generalFrictionContactModel&
Foam::solid4GeneralContactFvPatchVectorField::frictionModel
(
    const label targetI
) const
{
    if (!localSlave()[targetI])
    {
        FatalErrorIn("frictionModel(const label targetI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (frictionModels_.empty())
    {
        calcFrictionModels(dict_);
    }

    return frictionModels_[targetI];
}


void Foam::solid4GeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneIndex() const
{
    if (zoneIndex_ != -1)
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneIndex()"
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
            "void Foam::solid4GeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneIndex()"
            "const"
        )   << "Face zone name " << zoneName
            << " not found.  Please check your zone definition." << nl
            << "Current faceZones are:" << mesh.faceZones().names()
            << abort(FatalError);
    }

    zoneIndex_ = zone.index();
}


void Foam::solid4GeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZone() const
{
    if (globalPolyPatchFaceZonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZone() const"
        )   << "globalPolyPatchFaceZonePtr_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    globalPolyPatchFaceZonePtr_ =
        new globalPolyPatch
		(
        patch().name(),
        patch().boundaryMesh().mesh()
		);
		        
		/*
		(
            mesh.faceZones()[globalPolyPatchFaceZoneIndex()]().localFaces(),
            mesh.faceZones()[globalPolyPatchFaceZoneIndex()]().localPoints()
        );
		*/
}


Foam::label Foam::solid4GeneralContactFvPatchVectorField::globalPolyPatchFaceZoneIndex() const
{
    if (zoneIndex_ == -1)
    {
        calcGlobalPolyPatchFaceZoneIndex();
    }

    return zoneIndex_;
}


const Foam::globalPolyPatch&
Foam::solid4GeneralContactFvPatchVectorField::globalPolyPatchFaceZone() const
{
    if (!globalPolyPatchFaceZonePtr_)
    {
        calcGlobalPolyPatchFaceZone();
    }

    return *globalPolyPatchFaceZonePtr_;
}


Foam::globalPolyPatch& Foam::solid4GeneralContactFvPatchVectorField::globalPolyPatchFaceZone()
{
    if (!globalPolyPatchFaceZonePtr_)
    {
        calcGlobalPolyPatchFaceZone();
    }

    return *globalPolyPatchFaceZonePtr_;
}


const Foam::globalPolyPatch&
Foam::solid4GeneralContactFvPatchVectorField::targetGlobalPolyPatchFaceZone
(
    const label targetI
) const
{
    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    const solid4GeneralContactFvPatchVectorField& targetPatchField =
        refCast<const solid4GeneralContactFvPatchVectorField>
        (
            field.boundaryField()[targetPatchIndices()[targetI]]
        );

    return targetPatchField.globalPolyPatchFaceZone();
}


Foam::globalPolyPatch&
Foam::solid4GeneralContactFvPatchVectorField::targetGlobalPolyPatchFaceZone
(
    const label targetI
)
{
    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    // Const cast away the const-ness
    solid4GeneralContactFvPatchVectorField& targetPatchField =
        const_cast<solid4GeneralContactFvPatchVectorField&>
        (
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()[targetPatchIndices()[targetI]]
            )
        );

    return targetPatchField.globalPolyPatchFaceZone();
}


void Foam::solid4GeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneToZones() const
{
    // Create zone-to-zone interpolation
    if (!zoneToZones_.empty())
    {
        FatalErrorIn
        (
            "void solid4GeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneToZones()"
            "const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    zoneToZones_.setSize(targetPatchNames().size());

    const boolList& locSlave = localSlave();

    forAll(zoneToZones_, targetI)
    {
        // Only the local slave creates the interpolator
        if (locSlave[targetI])
        {
            zoneToZones_.set
                (
                    targetI,
                    new newGgiStandAlonePatchInterpolation
                    (
                        targetGlobalPolyPatchFaceZone(targetI).globalPatch(), // master
                        globalPolyPatchFaceZone().globalPatch(), // slave
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
        }
    }
}


const Foam::newGgiStandAlonePatchInterpolation&
Foam::solid4GeneralContactFvPatchVectorField::zoneToZone
(
    const label targetI
) const
{
    if (!localSlave()[targetI])
    {
        FatalErrorIn("zoneToZone(const label targetI)")
            << "Only the local slave can call the zoneToZone interpolator"
            << abort(FatalError);
    }

    if (zoneToZones_.empty())
    {
        word zoneName =
            patch().boundaryMesh().mesh().faceZones()[globalPolyPatchFaceZoneIndex()].name();

        if (debug)
        {
            Info<< "Initializing the GGI interpolators for " << zoneName
                << endl;
        }

        calcGlobalPolyPatchFaceZoneToZones();
    }

    return zoneToZones_[targetI];
}


Foam::newGgiStandAlonePatchInterpolation&
Foam::solid4GeneralContactFvPatchVectorField::zoneToZone(const label targetI)
{
    if (!localSlave()[targetI])
    {
        FatalErrorIn("zoneToZone(const label targetI)")
            << "Only the local slave can call the zoneToZone interpolator"
            << abort(FatalError);
    }

    if (zoneToZones_.empty())
    {
        word zoneName =
            patch().boundaryMesh().mesh().faceZones()[globalPolyPatchFaceZoneIndex()].name();

        if (debug)
        {
            Info<< "Initializing the GGI interpolators for " << zoneName
                << endl;
        }

        calcGlobalPolyPatchFaceZoneToZones();
    }

    return zoneToZones_[targetI];
}


Foam::label Foam::solid4GeneralContactFvPatchVectorField::findTargetID
(
    const label patchID
) const
{
    label targetI = -1;

    const labelList targetIDs = targetPatchIndices();

    forAll(targetIDs, I)
    {
        if (patchID == targetIDs[I])
        {
            targetI = I;
            break;
        }
    }

    if (targetI == -1)
    {
        FatalErrorIn("findTargetID(const label patchID)")
            << "target patch not found!" << abort(FatalError);
    }

    return targetI;
}


void Foam::solid4GeneralContactFvPatchVectorField::makeCurPatchTractions() const
{
    if (curPatchTractionPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::"
            "makeCurPatchTractions() const"
        )   << "curPatchTractionPtr_ already set" << abort(FatalError);
    }

    curPatchTractionPtr_ =
        new List<vectorField>
        (
            targetPatchNames().size(),
            vectorField(patch().size(), vector::zero)
        );
}


void Foam::solid4GeneralContactFvPatchVectorField::calcQc() const
{
    if (QcPtr_)
    {
        FatalErrorIn("solid4GeneralContactFvPatchVectorField::calcQc")
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


    // Accumulate Qc for from all targets

    const scalar deltaT = patch().boundaryMesh().mesh().time().deltaT().value();

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    const boolList& locSlave = localSlave();

    forAll(locSlave, targetI)
    {
        vectorField curPatchSlip(patch().size(), vector::zero);

        // Calculate slip
        if (locSlave[targetI])
        {
            curPatchSlip = frictionModel(targetI).slip();
        }
        else
        {
            const solid4GeneralContactFvPatchVectorField& targetPatchField =
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()[targetPatchIndices()[targetI]]
            );

            const label locTargetID =
                targetPatchField.findTargetID(patch().index());

            vectorField targetPatchSlip =
                targetPatchField.frictionModel(locTargetID).slip();

            vectorField targetGlobalPolyPatchFaceZoneSlip = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(targetPatchSlip);
			/*
                zoneField
                (
                    targetGlobalPolyPatchFaceZoneIndices()[targetI],
                    targetPatchIndices()[targetI],
                    targetPatchSlip
                );
				*/

            // Interpolate from target to the current patch
            // Face-to-face
			
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;

            vectorField curZoneSlip =
                targetPatchField.zoneToZone(locTargetID).slaveToMaster
                (
                    targetGlobalPolyPatchFaceZoneSlip
                );

            curPatchSlip = globalPolyPatchFaceZone().globalFaceToPatch(curZoneSlip);
			/*
                patchField
                (
                    patch().index(),
                    globalPolyPatchFaceZoneIndex(),
                    curZoneSlip
                );
				*/
        }

        // Heat flux rate: rate of dissipated frictional energy
        // The dot product of the traction vectors and the slip vectors
        // gives the dissipated frictional energy rate per unit area, which
        // is always positive
        Qc += mag(curTraction & (curPatchSlip/deltaT));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::solid4GeneralContactFvPatchVectorField::
solid4GeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localSlavePtr_(NULL),
    targetPatchNamesPtr_(NULL),
    targetPatchIndicesPtr_(NULL),
    zoneIndex_(-1),
    targetGlobalPolyPatchFaceZoneNamesPtr_(NULL),
    targetGlobalPolyPatchFaceZoneIndicesPtr_(NULL),
    rigidMaster_(false),
    dict_(NULL),
    normalModels_(0),
    frictionModels_(0),
    globalPolyPatchFaceZonePtr_(NULL),
    zoneToZones_(0),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(0.0)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solid4GeneralContactFvPatchVectorField::"
            "solid4GeneralContactFvPatchVectorField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        ) << endl;
    }
}


Foam::solid4GeneralContactFvPatchVectorField::
solid4GeneralContactFvPatchVectorField
(
    const solid4GeneralContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localSlavePtr_(NULL),
    targetPatchNamesPtr_(NULL),
    targetPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    targetGlobalPolyPatchFaceZoneNamesPtr_(NULL),
    targetGlobalPolyPatchFaceZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    globalPolyPatchFaceZonePtr_(NULL),
    zoneToZones_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solid4GeneralContactFvPatchVectorField::"
            "solid4GeneralContactFvPatchVectorField"
            "("
            "    const solid4GeneralContactFvPatchVectorField& ptf,"
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF,"
            "    const fvPatchFieldMapper& mapper"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.firstPatchInListPtr_)
    {
        firstPatchInListPtr_ = new bool(*ptf.firstPatchInListPtr_);
    }

    if (ptf.firstPatchInListIndexPtr_)
    {
        firstPatchInListIndexPtr_ = new label(*ptf.firstPatchInListIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.targetPatchNamesPtr_)
    {
        targetPatchNamesPtr_ = new wordList(*ptf.targetPatchNamesPtr_);
    }

    if (ptf.targetPatchIndicesPtr_)
    {
        targetPatchIndicesPtr_ = new labelList(*ptf.targetPatchIndicesPtr_);
    }

    if (ptf.targetGlobalPolyPatchFaceZoneNamesPtr_)
    {
        targetGlobalPolyPatchFaceZoneNamesPtr_ = new wordList(*ptf.targetGlobalPolyPatchFaceZoneNamesPtr_);
    }

    if (ptf.targetGlobalPolyPatchFaceZoneIndicesPtr_)
    {
        targetGlobalPolyPatchFaceZoneIndicesPtr_ = new labelList(*ptf.targetGlobalPolyPatchFaceZoneIndicesPtr_);
    }
	
	/*
	//Leave this for the moment since this constructor isn't called 
    if (ptf.globalPolyPatchFaceZonePtr_)
    {
        globalPolyPatchFaceZonePtr_ = new standAlonePatch(*ptf.globalPolyPatchFaceZonePtr_);
    }
	*/

    if (!ptf.zoneToZones_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solid4GeneralContactFvPatchVectorField::"
            "solid4GeneralContactFvPatchVectorField"
            "("
            "    const solid4GeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solid4GeneralContact: zoneToZone GGI interpolators not mapped"
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


Foam::solid4GeneralContactFvPatchVectorField::
solid4GeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localSlavePtr_(NULL),
    targetPatchNamesPtr_(NULL),
    targetPatchIndicesPtr_(NULL),
    zoneIndex_(-1),
    targetGlobalPolyPatchFaceZoneNamesPtr_(NULL),
    targetGlobalPolyPatchFaceZoneIndicesPtr_(NULL),
    rigidMaster_(dict.lookupOrDefault<Switch>("rigidMaster", false)),
    dict_(dict),
    normalModels_(0),
    frictionModels_(0),
    globalPolyPatchFaceZonePtr_(0),
    zoneToZones_(0),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(0.0)
{
    Info<< "Creating " << solid4GeneralContactFvPatchVectorField::typeName
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


Foam::solid4GeneralContactFvPatchVectorField::
solid4GeneralContactFvPatchVectorField
(
    const solid4GeneralContactFvPatchVectorField& ptf
)
:
    solidTractionFvPatchVectorField(ptf),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localSlavePtr_(NULL),
    targetPatchNamesPtr_(NULL),
    targetPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    targetGlobalPolyPatchFaceZoneNamesPtr_(NULL),
    targetGlobalPolyPatchFaceZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    globalPolyPatchFaceZonePtr_(NULL),
    zoneToZones_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solid4GeneralContactFvPatchVectorField::"
            "solid4GeneralContactFvPatchVectorField"
            "("
            "    const solid4GeneralContactFvPatchVectorField& ptf"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.firstPatchInListPtr_)
    {
        firstPatchInListPtr_ = new bool(*ptf.firstPatchInListPtr_);
    }

    if (ptf.firstPatchInListIndexPtr_)
    {
        firstPatchInListIndexPtr_ = new label(*ptf.firstPatchInListIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.targetPatchNamesPtr_)
    {
        targetPatchNamesPtr_ = new wordList(*ptf.targetPatchNamesPtr_);
    }

    if (ptf.targetPatchIndicesPtr_)
    {
        targetPatchIndicesPtr_ = new labelList(*ptf.targetPatchIndicesPtr_);
    }

    if (ptf.targetGlobalPolyPatchFaceZoneNamesPtr_)
    {
        targetGlobalPolyPatchFaceZoneNamesPtr_ = new wordList(*ptf.targetGlobalPolyPatchFaceZoneNamesPtr_);
    }

    if (ptf.targetGlobalPolyPatchFaceZoneIndicesPtr_)
    {
        targetGlobalPolyPatchFaceZoneIndicesPtr_ = new labelList(*ptf.targetGlobalPolyPatchFaceZoneIndicesPtr_);
    }

	/*
	//Leave this for the moment since this constructor isn't called 
    if (ptf.globalPolyPatchFaceZonePtr_)
    {
        globalPolyPatchFaceZonePtr_ = new standAlonePatch(*ptf.globalPolyPatchFaceZonePtr_);
    }
	*/

    if (!ptf.zoneToZones_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solid4GeneralContactFvPatchVectorField::"
            "solid4GeneralContactFvPatchVectorField"
            "("
            "    const solid4GeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solid4GeneralContact: zoneToZone GGI interpolators not mapped"
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


Foam::solid4GeneralContactFvPatchVectorField::
solid4GeneralContactFvPatchVectorField
(
    const solid4GeneralContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localSlavePtr_(NULL),
    targetPatchNamesPtr_(NULL),
    targetPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    targetGlobalPolyPatchFaceZoneNamesPtr_(NULL),
    targetGlobalPolyPatchFaceZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    globalPolyPatchFaceZonePtr_(NULL),
    zoneToZones_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solid4GeneralContactFvPatchVectorField::"
            "solid4GeneralContactFvPatchVectorField"
            "("
            "    const solid4GeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.firstPatchInListPtr_)
    {
        firstPatchInListPtr_ = new bool(*ptf.firstPatchInListPtr_);
    }

    if (ptf.firstPatchInListIndexPtr_)
    {
        firstPatchInListIndexPtr_ = new label(*ptf.firstPatchInListIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.targetPatchNamesPtr_)
    {
        targetPatchNamesPtr_ = new wordList(*ptf.targetPatchNamesPtr_);
    }

    if (ptf.targetPatchIndicesPtr_)
    {
        targetPatchIndicesPtr_ = new labelList(*ptf.targetPatchIndicesPtr_);
    }

    if (ptf.targetGlobalPolyPatchFaceZoneNamesPtr_)
    {
        targetGlobalPolyPatchFaceZoneNamesPtr_ = new wordList(*ptf.targetGlobalPolyPatchFaceZoneNamesPtr_);
    }

    if (ptf.targetGlobalPolyPatchFaceZoneIndicesPtr_)
    {
        targetGlobalPolyPatchFaceZoneIndicesPtr_ = new labelList(*ptf.targetGlobalPolyPatchFaceZoneIndicesPtr_);
    }

    /*
	//Be careful: Leave this for the moment since this constructor is called 
    if (ptf.globalPolyPatchFaceZonePtr_)
    {
        globalPolyPatchFaceZonePtr_ = new standAlonePatch(*ptf.globalPolyPatchFaceZonePtr_);
    }
	*/

    if (!ptf.zoneToZones_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solid4GeneralContactFvPatchVectorField::"
            "solid4GeneralContactFvPatchVectorField"
            "("
            "    const solid4GeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solid4GeneralContact: zoneToZone GGI interpolators not mapped"
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


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::solid4GeneralContactFvPatchVectorField::
~solid4GeneralContactFvPatchVectorField()
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solid4GeneralContactFvPatchVectorField::"
            "~solid4GeneralContactFvPatchVectorField()"
        ) << endl;
    }

    deleteDemandDrivenData(firstPatchInListPtr_);
    deleteDemandDrivenData(firstPatchInListIndexPtr_);
    deleteDemandDrivenData(localSlavePtr_);
    deleteDemandDrivenData(targetPatchNamesPtr_);
    deleteDemandDrivenData(targetPatchIndicesPtr_);
    deleteDemandDrivenData(targetGlobalPolyPatchFaceZoneNamesPtr_);
    deleteDemandDrivenData(targetGlobalPolyPatchFaceZoneIndicesPtr_);

    normalModels_.clear();
    frictionModels_.clear();

    deleteDemandDrivenData(globalPolyPatchFaceZonePtr_);

    zoneToZones_.clear();

    deleteDemandDrivenData(curPatchTractionPtr_);
    deleteDemandDrivenData(QcPtr_);
    deleteDemandDrivenData(QcsPtr_);
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //


bool Foam::solid4GeneralContactFvPatchVectorField::firstPatchInList() const
{
    if (!firstPatchInListPtr_)
    {
        calcFirstPatchInList();
    }

    return *firstPatchInListPtr_;
}


Foam::label Foam::solid4GeneralContactFvPatchVectorField::firstPatchInListIndex()
const
{
    if (!firstPatchInListIndexPtr_)
    {
        calcFirstPatchInListIndex();
    }

    return *firstPatchInListIndexPtr_;
}


const Foam::List<Foam::word>&
Foam::solid4GeneralContactFvPatchVectorField::targetPatchNames() const
{
    if (!targetPatchNamesPtr_)
    {
        calctargetPatchNames();
    }

    return *targetPatchNamesPtr_;
}


const Foam::List<Foam::label>&
Foam::solid4GeneralContactFvPatchVectorField::targetPatchIndices() const
{
    if (!targetPatchIndicesPtr_)
    {
        calctargetPatchNames();
    }

    return *targetPatchIndicesPtr_;
}


const Foam::List<Foam::word>&
Foam::solid4GeneralContactFvPatchVectorField::targetGlobalPolyPatchFaceZoneNames() const
{
    if (!targetGlobalPolyPatchFaceZoneNamesPtr_)
    {
        calcTargetGlobalPolyPatchFaceZoneNames();
    }

    return *targetGlobalPolyPatchFaceZoneNamesPtr_;
}


const Foam::List<Foam::label>&
Foam::solid4GeneralContactFvPatchVectorField::targetGlobalPolyPatchFaceZoneIndices() const
{
    if (!targetGlobalPolyPatchFaceZoneIndicesPtr_)
    {
        calcTargetGlobalPolyPatchFaceZoneNames();
    }

    return *targetGlobalPolyPatchFaceZoneIndicesPtr_;
}


// Map from self
void Foam::solid4GeneralContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    FatalErrorIn
    (
        "void Foam::solid4GeneralContactFvPatchVectorField::autoMap"
        "("
        "    const fvPatchFieldMapper& m"
        ")"
    )   << "member mapping not implemented" << endl;

    solidTractionFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::solid4GeneralContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    FatalErrorIn
    (
        "void Foam::solid4GeneralContactFvPatchVectorField::rmap"
        "("
        "    const fvPatchField<vector>& ptf,"
        "    const labelList& addr"
        ")"
    )   << "member mapping not implemented" << endl;

    solidTractionFvPatchVectorField::rmap(ptf, addr);
}


void Foam::solid4GeneralContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
	
	#if(ISDEBUG)
	Info<< "Check 1: patch().name() in updateCoeffs() "<<patch().name()<<endl;
	Info<< "patch().index() in updateCoeffs() "<<patch().index()<<endl;	
	Info<< "targetPatchNames().size() in updateCoeffs() "<<targetPatchNames().size()<<endl;	
	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;	
	#endif

    boolList activeContactPairs(targetPatchNames().size(), true);

    // if it is a new time step then reset iCorr
    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Delete friction heat rate to force its recalculation when thermal
        // boundaries ask for it
        deleteDemandDrivenData(QcPtr_);
        deleteDemandDrivenData(QcsPtr_);

        if (firstPatchInList())
        {
            forAll(activeContactPairs, activePatchI)
            {
                // Let the contact models know that it is a new time-step, in
                // case they need to update anything
                normalModel(activePatchI).newTimeStep();
                frictionModel(activePatchI).newTimeStep();
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
        if (firstPatchInList())
        {
			#if(ISDEBUG)
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
			#endif
           // Move the master and slave zone to the deformed configuration
            moveZonesToDeformedConfiguration();
        }


        // Clear interpolator weights

        forAll(activeContactPairs, slaveI)
        {
            if (localSlave()[slaveI])
            {
                zoneToZone(slaveI).movePoints
                (
                    tensorField(0), tensorField(0), vectorField(0)
                );
            }
        }


        // Accumulated traction for the current patch
        vectorField curPatchTraction(patch().size(), vector::zero);

        // Only the local masters calculates the contact force and the local
        // master interpolates this force
        const boolList& locSlave = localSlave();

        // Create master bounding box used for quick check
        boundBox masterBb(globalPolyPatchFaceZone().patch().localPoints(), false);

        // The BB may have zero thickness in one of the directions e.g. for a
        // flat patch, so we will check for this and, if found, create an offset
        const scalar bbOff = bbOffset();
        if (masterBb.minDim() < bbOff)
        {
            const vector bbDiag = masterBb.max() - masterBb.min();

            if (bbDiag.x() < bbOff)
            {
                vector offset(bbOff, 0, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
            //else 
			if (bbDiag.y() < bbOff)
            {
                vector offset(0, bbOff, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
            //else 
			if (bbDiag.z() < bbOff)
            {
                vector offset(0, 0, bbOff);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
        }

        forAll(activeContactPairs, activePatchI)
        {
			#if(ISDEBUG)
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
			Info<<"activePatchI pair in updateCoeffs(): "<<activePatchI<<endl;
			#endif
			
            // Perform quick check to find potential contacting pairs
            // The quick check is based on the bounding box (BB) of the contact
            // pairs: if the BBs of pair intersect then we will designate the
            // pair as active.
			
			//***************** Start boundBox comment ************
			#if(ISDEBUG)
            // Create shadow bounding box
            boundBox shadowBb(targetGlobalPolyPatchFaceZone(activePatchI).patch().localPoints(), false);

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
                //else 
				if (bbDiag.y() < bbOff)
                {
                    vector offset(0, bbOff, 0);
                    shadowBb.min() -= offset;
                    shadowBb.max() += offset;
                }
                //else 
				if (bbDiag.z() < bbOff)
                {
                    vector offset(0, 0, bbOff);
                    shadowBb.min() -= offset;
                    shadowBb.max() += offset;
                }
            }

            if (masterBb.overlaps(shadowBb))
            {
                activeContactPairs[activePatchI] = true;
            }
			
			//***************** End boundBox comment ************
			#endif

            // Call normal and frction contact models for active contacting
            // pairs
            // Accumulate contact force contributions for all active contact
            // pairs

            if (activeContactPairs[activePatchI])
            {
				#if(!ActivePairDEBUG)				
				Info<<"Checking mMASTER or sSLAVE in updateCoeffs() line:"<<__LINE__<<endl;
				Info<<"locSlave[activePatchI] in updateCoeffs(): "<<locSlave[activePatchI]<<endl;				
				Info<< "patch().name() in updateCoeffs() "<<patch().name()<<endl;	
				Info<< "patch INDEX in updateCoeffs() "<<patch().index()<<endl;
				Info<< "shadowPatch INDEX in updateCoeffs() "<<targetPatchIndices()[activePatchI]<<endl;
				Info<<"shadowPatchNames in updateCoeffs() "<<targetPatchNames()[activePatchI]<<endl;
				#endif
                if (locSlave[activePatchI])
                {
					#if(ActivePairDEBUG)
					Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
					Info<<"globalPolyPatchFaceZone().globalPatch() in updateCoeffs(): "<<globalPolyPatchFaceZone().globalPatch()<<endl;
                    Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
					#endif
					// Correct normal and friction contact models for the
                    // current contact pair

                    // Calculate the slave patch face unit normals as they are
                    // units by both the normal and friction models
                    const vectorField patchFaceNormals = 
										globalPolyPatchFaceZone().globalFaceToPatch
										(
										globalPolyPatchFaceZone().globalPatch().faceNormals()
										);
					#if(ISDEBUG)
					Info<<"slave patchFaceNormals in updateCoeffs(): "<<patchFaceNormals<<endl;
					#endif
					 /*
                        patchField
                        (
                            targetPatchIndices()[activePatchI],
                            shadowZoneIndices()[activePatchI],
                            targetGlobalPolyPatchFaceZone(activePatchI).globalPatch().faceNormals()
                        );
						*/
						
						#if(ISDEBUG)
						Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
						#endif
                    // Interpolate the master displacement increment to the
                    // slave patch as it is required by specific normal and
                    // friction contact models

                    vectorField patchDD(patch().size(), vector::zero);
                    vectorField targetPatchDD
                    (
                        patch().boundaryMesh()
                        [
                            targetPatchIndices()[activePatchI]
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
                        targetPatchDD =
                            DD.boundaryField()[targetPatchIndices()[activePatchI]];
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
                        targetPatchDD =
                            D.boundaryField()[targetPatchIndices()[activePatchI]]
                          - D.oldTime().boundaryField()
                            [
                                targetPatchIndices()[activePatchI]
                            ];
                    }
					
					#if(ISDEBUG)
					Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
					Info<<"patchDD in updateCoeffs(): "<<patchDD<<endl;
					Info<<"targetPatchDD in updateCoeffs() line:"<<targetPatchDD<<endl;
					#endif
					
                    // Master zone DD
                    const vectorField tagetZoneDD = targetGlobalPolyPatchFaceZone(activePatchI).patchFaceToGlobal(targetPatchDD);
					
					#if(ISDEBUG)
					Info<<"tagetZoneDD in updateCoeffs(): "<<tagetZoneDD<<endl;
					#endif
					
					/*	
					   zoneField
                        (
                            zoneIndex(),
                            patch().index(),
                            patchDD
                        );
					*/
					
					#if(ISDEBUG)
					Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
					#endif
					
                    // Master patch DD interpolated to the slave patch
                    const vectorField patchDDInterpToShadowPatch =
						globalPolyPatchFaceZone().globalFaceToPatch
										(
										zoneToZone(activePatchI).masterToSlave(tagetZoneDD)()
										);
					
					#if(ISDEBUG)
					Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
					
					Info<<"patchDDInterpToShadowPatch in updateCoeffs(): "<<patchDDInterpToShadowPatch<<endl;
					#endif
					/*
                        patchField
                        (
                            targetPatchIndices()[activePatchI],
                            shadowZoneIndices()[activePatchI],
                            zoneToZone(activePatchI).masterToSlave(zoneDD)()
                        );
						*/

                    //FatalError
                    //    << "Disabled: use jasakSolidContact" << abort(FatalError);
                     
					 normalModel(activePatchI).correct
                     (
                         patchFaceNormals,
                         //zoneToZone(activePatchI),
						 //targetGlobalPolyPatchFaceZone(activePatchI).globalPointToPatch
						 //(
						 zoneToZone(activePatchI).slavePointDistanceToIntersection(),
						 //),
                         targetPatchDD,
                         patchDDInterpToShadowPatch
                     );
					 
					 #if(ISDEBUG)
					 Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
					 Info<<"normalModel(activePatchI).slavePressure()in updateCoeffs(): "<<normalModel(activePatchI).slavePressure()<<endl;
					 #endif
					 
                    frictionModel(activePatchI).correct
                    (
                        normalModel(activePatchI).slavePressure(),
                        patchFaceNormals,
                        normalModel(activePatchI).areaInContact(),
                        targetPatchDD,
                        patchDDInterpToShadowPatch
                    );
					
					#if(ISDEBUG)
					Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
					#endif
					
                    // Accumulate traction

                    curPatchTractions(activePatchI) =
                        frictionModel(activePatchI).slaveTraction()
                        + normalModel(activePatchI).slavePressure();

                    curPatchTraction += curPatchTractions(activePatchI);
					
					#if(!ISDEBUG)
					Info<<"SLAVE curPatchTraction in updateCoeffs() :"<<curPatchTraction<<endl;
					#endif
				}
                else // local master
                {
					#if(ActivePairDEBUG)
					Info<<"MASTER in updateCoeffs() line:"<<__LINE__<<endl;
                    #endif
					
					// Get traction from local slave

                    const volVectorField& field =
                        db().lookupObject<volVectorField>
                        (
                            dimensionedInternalField().name()
                        );

                    const solid4GeneralContactFvPatchVectorField&
                        localMasterField =
                        refCast<const solid4GeneralContactFvPatchVectorField>
                        (
                            field.boundaryField()
                            [
                                targetPatchIndices()[activePatchI]
                            ]
                        );

                    const label masterTargetI =
                        localMasterField.findTargetID(patch().index());

                    vectorField shadowPatchTraction =
                        -localMasterField.frictionModel
                        (
                            masterTargetI
                        ).slaveTractionForMaster()
                        -localMasterField.normalModel
                        (
                            masterTargetI
                        ).slavePressure();

                    vectorField shadowZoneTraction = targetGlobalPolyPatchFaceZone(activePatchI).patchFaceToGlobal(shadowPatchTraction);
                        /*
						zoneField
                        (
                            shadowZoneIndices()[activePatchI],
                            targetPatchIndices()[activePatchI],
                            shadowPatchTraction
                        );
						*/
					
					#if(ISDEBUG)
					Info<<"MASTER in updateCoeffs() line:"<<__LINE__<<endl;
					#endif
					
                    // Face-to-face
                    vectorField masterZoneTraction =
                        localMasterField.zoneToZone
                        (
                            masterTargetI
                        ).slaveToMaster(shadowZoneTraction);

                    // We store master patch traction as thermalGeneralContact
                    // uses it
                    curPatchTractions(activePatchI) =
						globalPolyPatchFaceZone().globalFaceToPatch(masterZoneTraction);
						/*
                        patchField
                        (
                            patch().index(),
                            zoneIndex(),
                            masterZoneTraction
                        );
						*/

                    curPatchTraction += curPatchTractions(activePatchI);
					
					#if(!ISDEBUG)
					Info<<"MASTER curPatchTraction in updateCoeffs() :"<<curPatchTraction<<endl;
					#endif
                }
            } // if contact pair is active
        } // forAll contact pairs

        // Set master gradient based on accumulated traction
        traction() = curPatchTraction;
    }

    solidTractionFvPatchVectorField::updateCoeffs();
	
	Info<<"In updateCoeffs() line"<<__LINE__<<endl;
		//******************** Testing face coordinates ************************
		const fvMesh& mesh = patch().boundaryMesh().mesh();
		label PatchI = mesh.boundaryMesh().findPatchID("patchName");
		// const vectorField& faceCentres = mesh.boundaryMesh()[patchI].faceCentres();
		Info<<"In updateCoeffs() line"<<__LINE__<<endl;
		
		const fvPatchVectorField& faceCentres = mesh.C().boundaryField()[PatchI];
		
		const vectorField& patchFaceCenter = patch().Cf();
		Info<<"In updateCoeffs() line"<<__LINE__<<endl;

        //const scalarField& zpatchFaceCenter = patchFaceCenter.component(2);	
		
        Info<< "Field patchFaceCenter: " << patchFaceCenter << endl;
        //Info<< "Field zpatchFaceCenter:" << zpatchFaceCenter << endl;		

        //label celli(0);
        //const scalar z0 = patchFaceCenter[celli].z();
        //Info<< "Just the first cell: " << z0 << endl;

        Info<< "All cells Cf():" << endl;
        forAll(patchFaceCenter, facei) 
		{
            Info<< "patchFaceCenter ["<<facei<<"]:" << patchFaceCenter[facei] << endl;
        }
		//******************** End testing face coordinates *********************
		Info<<"In updateCoeffs() line"<<__LINE__<<endl;	
}

/*
void Foam::solid4GeneralContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    boolList activeContactPairs(targetPatchNames().size(), false);

    // if it is a new time step then reset iCorr
    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Delete friction heat rate to force its recalculation when thermal
        // boundaries ask for it
        deleteDemandDrivenData(QcPtr_);
        deleteDemandDrivenData(QcsPtr_);

        if (firstPatchInList())
        {
            forAll(activeContactPairs, targetI)
            {
                // Let the contact models know that it is a new time-step, in
                // case they need to update anything
                normalModel(targetI).newTimeStep();
                frictionModel(targetI).newTimeStep();
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
        if (firstPatchInList())
        {
           // Move the master and slave zone to the deformed configuration
            moveZonesToDeformedConfiguration();
        }


        // Clear interpolator weights

        forAll(activeContactPairs, slaveI)
        {
            if (localSlave()[slaveI])
            {
                zoneToZone(slaveI).movePoints
                (
                    tensorField(0), tensorField(0), vectorField(0)
                );
            }
        }


        // Accumulated traction for the current patch
        vectorField curPatchTraction(patch().size(), vector::zero);

        // Only the local masters calculates the contact force and the local
        // master interpolates this force
        const boolList& locSlave = localSlave();

        // Create master bounding box used for quick check
        boundBox masterBb(globalPolyPatchFaceZone().patch().localPoints(), false);

        // The BB may have zero thickness in one of the directions e.g. for a
        // flat patch, so we will check for this and, if found, create an offset
        const scalar bbOff = bbOffset();
        if (masterBb.minDim() < bbOff)
        {
            const vector bbDiag = masterBb.max() - masterBb.min();

            if (bbDiag.x() < bbOff)
            {
                vector offset(bbOff, 0, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
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

        forAll(activeContactPairs, targetI)
        {
            // Perform quick check to find potential contacting pairs
            // The quick check is based on the bounding box (BB) of the contact
            // pairs: if the BBs of pair intersect then we will designate the
            // pair as active.

            // Create target bounding box
            boundBox targetBb(targetGlobalPolyPatchFaceZone(targetI).patch().localPoints(), false);

            // Check for a zero dimension in the targetBb
            if (targetBb.minDim() < bbOff)
            {
                const vector bbDiag = targetBb.max() - targetBb.min();

                if (bbDiag.x() < bbOff)
                {
                    vector offset(bbOff, 0, 0);
                    targetBb.min() -= offset;
                    targetBb.max() += offset;
                }
                else if (bbDiag.y() < bbOff)
                {
                    vector offset(0, bbOff, 0);
                    targetBb.min() -= offset;
                    targetBb.max() += offset;
                }
                else if (bbDiag.z() < bbOff)
                {
                    vector offset(0, 0, bbOff);
                    targetBb.min() -= offset;
                    targetBb.max() += offset;
                }
            }

            if (masterBb.overlaps(targetBb))
            {
                activeContactPairs[targetI] = true;
            }

            // Call normal and frction contact models for active contacting
            // pairs
            // Accumulate contact force contributions for all active contact
            // pairs

            if (activeContactPairs[targetI])
            {
                if (locSlave[targetI])
                {
                    // Correct normal and friction contact models for the
                    // current contact pair

                    // Calculate the slave patch face unit normals as they are
                    // units by both the normal and friction models
                    const vectorField targetPatchFaceNormals = 
										targetGlobalPolyPatchFaceZone(targetI).globalFaceToPatch
										(
										targetGlobalPolyPatchFaceZone(targetI).globalPatch().faceNormals()
										);

                    // Interpolate the master displacement increment to the
                    // slave patch as it is required by specific normal and
                    // friction contact models

                    vectorField patchDD(patch().size(), vector::zero);
                    vectorField targetPatchDD
                    (
                        patch().boundaryMesh()
                        [
                            targetPatchIndices()[targetI]
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
                        targetPatchDD =
                            DD.boundaryField()[targetPatchIndices()[targetI]];
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
                        targetPatchDD =
                            D.boundaryField()[targetPatchIndices()[targetI]]
                          - D.oldTime().boundaryField()
                            [
                                targetPatchIndices()[targetI]
                            ];
                    }

                    // Master zone DD
                    const vectorField zoneDD = globalPolyPatchFaceZone().patchFaceToGlobal(patchDD);
										
					Info<<"In updateCoeffs() line:"<<__LINE__<<endl;

                    // Master patch DD interpolated to the slave patch
                    const vectorField patchDDInterpToTargetPatch =
						targetGlobalPolyPatchFaceZone(targetI).globalFaceToPatch
										(
										zoneToZone(targetI).slaveToMaster(zoneDD)()
										);
						
                    
					 normalModel(targetI).correct
                     (
                         targetPatchFaceNormals,
                         //zoneToZone(targetI),
						 //targetGlobalPolyPatchFaceZone(targetI).globalPointToPatch
						 //(
						 zoneToZone(targetI).slavePointDistanceToIntersection(),
						 //),
                         targetPatchDD,
                         patchDDInterpToTargetPatch
                     );
					 
					 Info<<"In updateCoeffs() line:"<<__LINE__<<endl;

                    frictionModel(targetI).correct
                    (
                        normalModel(targetI).slavePressure(),
                        targetPatchFaceNormals,
                        normalModel(targetI).areaInContact(),
                        targetPatchDD,
                        patchDDInterpToTargetPatch
                    );
					
					Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
					
                    // Accumulate traction

                    curPatchTractions(targetI) =
                        frictionModel(targetI).slaveTraction()
                        + normalModel(targetI).slavePressure();

                    curPatchTraction += curPatchTractions(targetI);
                }
                else // local master
                {
                    // Get traction from local slave

                    const volVectorField& field =
                        db().lookupObject<volVectorField>
                        (
                            dimensionedInternalField().name()
                        );

                    const solid4GeneralContactFvPatchVectorField&
                        localMasterField =
                        refCast<const solid4GeneralContactFvPatchVectorField>
                        (
                            field.boundaryField()
                            [
                                targetPatchIndices()[targetI]
                            ]
                        );

                    const label masterTargetI =
                        localMasterField.findTargetID(patch().index());

                    vectorField targetPatchTraction =
                        -localMasterField.frictionModel
                        (
                            masterTargetI
                        ).slaveTractionForMaster()
                        -localMasterField.normalModel
                        (
                            masterTargetI
                        ).slavePressure();

                    vectorField targetGlobalPolyPatchFaceZoneTraction = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(targetPatchTraction);
                        					
					Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
					
                    // Face-to-face
                    vectorField masterZoneTraction =
                        localMasterField.zoneToZone
                        (
                            masterTargetI
                        ).slaveToMaster(targetGlobalPolyPatchFaceZoneTraction);

                    // We store master patch traction as thermalGeneralContact
                    // uses it
                    curPatchTractions(targetI) =
						globalPolyPatchFaceZone().globalFaceToPatch(masterZoneTraction);
						
                    curPatchTraction += curPatchTractions(targetI);
                }
            } // if contact pair is active
        } // forAll contact pairs

        // Set master gradient based on accumulated traction
        traction() = curPatchTraction;
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}
*/

const Foam::scalarField& Foam::solid4GeneralContactFvPatchVectorField::Qc() const
{
    if (!QcPtr_)
    {
        calcQc();
    }

    return *QcPtr_;
}


//- Increment of dissipated energy due to friction for each pair
const Foam::scalarField& Foam::solid4GeneralContactFvPatchVectorField::Qc
(
    const label targetI
) const
{
    if (!QcsPtr_)
    {
        calcQcs();
    }

    return (*QcsPtr_)[targetI];
}


void Foam::solid4GeneralContactFvPatchVectorField::calcQcs() const
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

    forAll(locSlave, targetI)
    {
        scalarField& Qc = (*QcsPtr_)[targetI];
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

        // Calculate Qc for targetI

        vectorField curPatchSlip(Qc.size(), vector::zero);

        // Calculate slip
        if (locSlave[targetI])
        {
            curPatchSlip = frictionModel(targetI).slip();
        }
        else
        {
            const solid4GeneralContactFvPatchVectorField& targetPatchField =
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()[targetPatchIndices()[targetI]]
            );

            const label locTargetID =
            targetPatchField.findTargetID(patch().index());

            vectorField targetPatchSlip =
            targetPatchField.frictionModel(locTargetID).slip();

            vectorField targetGlobalPolyPatchFaceZoneSlip = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(targetPatchSlip);
            /*
			zoneField
            (
                targetGlobalPolyPatchFaceZoneIndices()[targetI],
                targetPatchIndices()[targetI],
                targetPatchSlip
            );
			*/

            // Interpolate from target to the current patch
            // Face-to-face
			
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
			
            vectorField curZoneSlip =
            targetPatchField.zoneToZone(locTargetID).slaveToMaster
            (
                targetGlobalPolyPatchFaceZoneSlip
            );

            curPatchSlip = globalPolyPatchFaceZone().globalFaceToPatch(curZoneSlip);
            /*
			patchField
            (
                patch().index(),
                globalPolyPatchFaceZoneIndex(),
                curZoneSlip
            );
			*/
        }

        // Heat flux rate: rate of dissipated frictional energy
        // The dot product of the traction vectors and the slip vectors
        // gives the dissipated frictional energy rate per unit area, which
        // is always positive
        Qc = mag(curTraction & (curPatchSlip/deltaT));
    }
}


void Foam::solid4GeneralContactFvPatchVectorField::calcBbOffset() const
{
    if (bbOffset_ != 0)
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::solid4GeneralContactFvPatchVectorField::"
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



Foam::scalar Foam::solid4GeneralContactFvPatchVectorField::bbOffset() const
{
    if (bbOffset_ == 0)
    {
        calcBbOffset();
    }

    return bbOffset_;
}


const Foam::vectorField&
Foam::solid4GeneralContactFvPatchVectorField::curPatchTractions
(
    const label targetI
) const
{
    if (!curPatchTractionPtr_)
    {
        makeCurPatchTractions();
    }

    return (*curPatchTractionPtr_)[targetI];
}


Foam::vectorField&
Foam::solid4GeneralContactFvPatchVectorField::curPatchTractions
(
    const label targetI
)
{
    if (!curPatchTractionPtr_)
    {
        makeCurPatchTractions();
    }

    return (*curPatchTractionPtr_)[targetI];
}



Foam::scalarField
Foam::solid4GeneralContactFvPatchVectorField::scaleTractionField() const
{
    if (scaleTractionFieldPtr_.empty())
    {
        makeScaleTractionField();
    }

    return scaleTractionFieldPtr_();
}


void Foam::solid4GeneralContactFvPatchVectorField::makeScaleTractionField() const
{
    if (scaleTractionFieldPtr_.valid())
    {
        FatalErrorIn(type() + "::makeScaleTractionField()")
            << "Pointer already set!" << abort(FatalError);
    }

    scaleTractionFieldPtr_.set(new scalarField(patch().size(), 1.0));
    scalarField& scaleTractionField = scaleTractionFieldPtr_();

    // Find all faces on the patch that are adjacent to faces on the
    // downstream patch

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const scalar scaleFactor =
        readScalar(dict_.lookup("downstreamScaleFactor"));

    // Downstream patch name
    const word patchName = dict_.lookup("downstreamPatchName");
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    if (patchID == -1)
    {
        FatalErrorIn(type() + "::makeScaleTractionField()")
            << "Cannot find patch " << patchName << abort(FatalError);
    }

    const unallocLabelList& faceCells = patch().faceCells();
    const cellList& cells = mesh.cells();

    forAll(scaleTractionField, fI)
    {
        const label cellID = faceCells[fI];
        const cell& curCell = cells[cellID];

        forAll(curCell, cfI)
        {
            const label cellFaceID = curCell[cfI];

            if (!mesh.isInternalFace(cellFaceID))
            {
                if (mesh.boundaryMesh().whichPatch(cellFaceID) == patchID)
                {
                    scaleTractionField[fI] = scaleFactor;
                    break;
                }
            }
        }
    }
}


Foam::PtrList<Foam::scalarField>&
Foam::solid4GeneralContactFvPatchVectorField::contactPerSlave()
{
    if (contactPerSlave_.size() == 0)
    {
        calcContactPerSlave();
    }

    return contactPerSlave_;
}


const Foam::PtrList<Foam::scalarField>&
Foam::solid4GeneralContactFvPatchVectorField::contactPerSlave() const
{
    if (contactPerSlave_.size() == 0)
    {
        calcContactPerSlave();
    }

    return contactPerSlave_;
}

void Foam::solid4GeneralContactFvPatchVectorField::calcContactPerSlave() const
{
    if (contactPerSlave_.size() > 0)
    {
        FatalErrorIn
        (
            "void thermalContactFvPatchScalarField::"
            "calcContactPerSlave() const"
        )   << "already calculated"
            << abort(FatalError);
    }

    contactPerSlave_.setSize(targetPatchNames().size());

    forAll(contactPerSlave_, i)
    {
        contactPerSlave_.set
        (
            i,
            new scalarField(patch().size(), 0.0)
        );
    }
}



// Write
void Foam::solid4GeneralContactFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("rigidMaster")
        << rigidMaster_ << token::END_STATEMENT << nl;

    // Write the dict from the first contact model

    const label targetI = 0;

    if (localSlave()[targetI])
    {
        os.writeKeyword("normalContactModel")
            << normalModel(targetI).type() << token::END_STATEMENT << nl;
        normalModel(targetI).writeDict(os);

        os.writeKeyword("frictionContactModel")
            << frictionModel(targetI).type() << token::END_STATEMENT << nl;
        frictionModel(targetI).writeDict(os);
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                dimensionedInternalField().name()
            );

        const solid4GeneralContactFvPatchVectorField& localSlaveField =
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()
                [
                    targetPatchIndices()[targetI]
                ]
            );

        const label localSlaveID =
            localSlaveField.findTargetID(patch().index());

        os.writeKeyword("normalContactModel")
            << localSlaveField.normalModel(localSlaveID).type()
            << token::END_STATEMENT << nl;
        localSlaveField.normalModel(localSlaveID).writeDict(os);

        os.writeKeyword("frictionContactModel")
            << localSlaveField.frictionModel(localSlaveID).type()
            << token::END_STATEMENT << nl;
        localSlaveField.frictionModel(localSlaveID).writeDict(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField, solid4GeneralContactFvPatchVectorField
    )
}


// ************************************************************************* //
