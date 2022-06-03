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
    solidGeneralContactFvPatchVectorField

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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


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

    // Shadow patch and zone indices
    const labelList& shadPatchIndices = targetPatchIndices();
    const labelList& shadZoneIndices = targetGlobalPolyPatchFaceZoneIndices();

    forAll(shadPatchIndices, targetI)
    {
        // Assemble the zone face displacement field to move the zones
        vectorField zoneD(globalPolyPatchFaceZone().globalPatch().size(), vector::zero);
        vectorField shadowZoneD(targetGlobalPolyPatchFaceZone(targetI).globalPatch().size(), vector::zero);

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
                DD.boundaryField()[shadPatchIndices[targetI]];

            zoneD = globalPolyPatchFaceZone().patchFaceToGlobal(patchDD);		
            //    zoneField(globalPolyPatchFaceZoneIndex(), patch().index(), patchDD);
			
            shadowZoneD = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(targetPatchDD);
			
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
                D.boundaryField()[shadPatchIndices[targetI]];

            zoneD = globalPolyPatchFaceZone().patchFaceToGlobal(patchD);
                //zoneField(globalPolyPatchFaceZoneIndex(), patch().index(), patchD);
            shadowZoneD = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(targetPatchDD);
                
        }

        // Interpolate the zone face field to the zone points
        const pointField zonePointD = globalPolyPatchFaceZone().interpolator().faceToPointInterpolate(zoneD);
            //zoneFaceToPointInterpolate(globalPolyPatchFaceZoneIndex(), zoneD, -1);
        const pointField shadowZonePointD =
			targetGlobalPolyPatchFaceZone(targetI).interpolator().faceToPointInterpolate
            (
              shadowZoneD
            );
			

        // The zone deformed points are the initial position plus the
        // displacement
        const pointField zoneNewPoints =
            mesh.faceZones()[globalPolyPatchFaceZoneIndex()]().localPoints()
          + zonePointD;
        const pointField shadowZoneNewPoints =
            mesh.faceZones()[shadZoneIndices[targetI]]().localPoints()
          + shadowZonePointD;

        // Move the zones

        // Remove zones weights
        if (targetI == 0)
        {
            globalPolyPatchFaceZone().movePoints(zoneNewPoints);
        }
        targetGlobalPolyPatchFaceZone(targetI).movePoints(shadowZoneNewPoints);

        // We need to use const_cast to move the standAlonePatch points as the
        // movePoints function only clears weights
        // Also, be careful to move the points are opposed to the localPoints
        if (targetI == 0)
        {
            const_cast<pointField&>(globalPolyPatchFaceZone().globalPatch().points()) = zoneNewPoints;
        }
        const_cast<pointField&>(targetGlobalPolyPatchFaceZone(targetI).globalPatch().points()) =
            shadowZoneNewPoints;
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcFirstPatchInList() const
{
    if (firstPatchInListPtr_)
    {
        FatalErrorIn
            (
                "void Foam::solidGeneralContactFvPatchVectorField::"
                "calcFirstPatchInList() const"
            )   << "firstPatchInListPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solidGeneralContact patch i.e. the one
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


void Foam::solidGeneralContactFvPatchVectorField::calcFirstPatchInListIndex() const
{
    if (firstPatchInListIndexPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcFirstPatchInListIndex() const"
        )   << "firstPatchInListIndexPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solidGeneralContact patch i.e. the one
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
            "calcFirstPatchInListIndex() const"
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
    if (localMasterPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcLocalSlave() const"
        )   << "localMasterPtr_ already set" << abort(FatalError);
    }

    localMasterPtr_ = new boolList(targetPatchNames().size(), false);

    boolList& localMaster = *localMasterPtr_;

    forAll(localMaster, targetI)
    {
        if (patch().index() < targetPatchIndices()[targetI])
        {
            localMaster[targetI] = true;

            Info<< "solidGeneralContact: "
                << targetPatchNames()[targetI] << " (master)" << " to "
                << patch().name() << " (slave)" << endl;
        }
    }
}


const Foam::boolList&
Foam::solidGeneralContactFvPatchVectorField::localMaster() const
{
    if (!localMasterPtr_)
    {
        calcLocalSlave();
    }

    return *localMasterPtr_;
}


void Foam::solidGeneralContactFvPatchVectorField::calcTargetPatchNames() const
{
    if (targetPatchNamesPtr_ || targetPatchIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcTargetPatchNames() const"
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

    targetPatchNamesPtr_ = new wordList(nShadPatches);
    wordList& shadowPatchNames = *targetPatchNamesPtr_;

    targetPatchIndicesPtr_ = new labelList(nShadPatches);
    labelList& shadowPatchIndices = *targetPatchIndicesPtr_;

    // Record shadow patch names

    label targetI = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solidGeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
            shadowPatchNames[targetI] = patch().boundaryMesh()[patchI].name();

            shadowPatchIndices[targetI++] = patchI;
        }
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcTargetGlobalPolyPatchFaceZoneNames() const
{
    if (targetGlobalPolyPatchFaceZoneNamesPtr_ || targetGlobalPolyPatchFaceZoneIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcTargetGlobalPolyPatchFaceZoneNames() const"
        )   << "shadowZoneNames_ or shadowZoneIndices_ already set"
            << abort(FatalError);
    }

    const wordList& shadNames = targetPatchNames();

    targetGlobalPolyPatchFaceZoneNamesPtr_ = new wordList(shadNames.size());
    wordList& shadowZoneNames = *targetGlobalPolyPatchFaceZoneNamesPtr_;

    targetGlobalPolyPatchFaceZoneIndicesPtr_ = new labelList(shadNames.size());
    labelList& shadowZoneIndices = *targetGlobalPolyPatchFaceZoneIndicesPtr_;

    const fvMesh& mesh = patch().boundaryMesh().mesh();
	

    forAll(shadNames, targetI)
    {
        word zoneName = shadNames[targetI] + "FaceZone";

        faceZoneID zone(zoneName, mesh.faceZones());

        if (!zone.active())
        {
            FatalErrorIn("solidGeneralContactFvPatchVectorField")
                << "Face zone name " << zoneName
                << " not found.  Please check your zone definition."
                << abort(FatalError);
        }

        shadowZoneNames[targetI] = zoneName;

        shadowZoneIndices[targetI] = zone.index();
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcNormalModels(const dictionary& dict) const
{
    if (!normalModels_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "calcNormalModel() const"
        )   << "normalModels already set" << abort(FatalError);
    }

    normalModels_.setSize(targetPatchNames().size());

    const boolList& locMaster = localMaster();
	

    forAll(normalModels_, targetI)
    {
		// Only the local slave creates the contact model
        if (locMaster[targetI])
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
                    patch(), //.boundaryMesh()[targetPatchIndices()[targetI]],
                    contactDict,
                    //targetPatchIndices()[targetI], // master
                    patch().index(), // master
					targetPatchIndices()[targetI], // slave
                    //targetGlobalPolyPatchFaceZone(targetI).globalPatch(), // master
                    globalPolyPatchFaceZone().globalPatch(), // master
					targetGlobalPolyPatchFaceZone(targetI).globalPatch() // slave
                )
            );
        }
    }
}


Foam::generalNormalContactModel&
Foam::solidGeneralContactFvPatchVectorField::normalModel(const label targetI)
{
    if (!localMaster()[targetI])
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
Foam::solidGeneralContactFvPatchVectorField::normalModel
(
    const label targetI
) const
{
    if (!localMaster()[targetI])
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


void Foam::solidGeneralContactFvPatchVectorField::calcFrictionModels(const dictionary& dict) const
{
    if (!frictionModels_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "calcFrictionModel(targetI) const"
        )   << "frictionModelPtr_[targetI] already set" << abort(FatalError);
    }

    frictionModels_.setSize(targetPatchNames().size());

    const boolList& locMaster = localMaster();

    forAll(frictionModels_, targetI)
    {
		if (locMaster[targetI])
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
                        patch(), //.boundaryMesh()[targetPatchIndices()[targetI]],
                        contactDict,
                        //targetPatchIndices()[targetI], // master
                        patch().index(), // slave
						targetPatchIndices()[targetI] // master
                    )
                );
        }
    }
}


Foam::generalFrictionContactModel&
Foam::solidGeneralContactFvPatchVectorField::frictionModel(const label targetI)
{
    if (!localMaster()[targetI])
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
Foam::solidGeneralContactFvPatchVectorField::frictionModel
(
    const label targetI
) const
{
    if (!localMaster()[targetI])
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


void Foam::solidGeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneIndex() const
{
    if (zoneIndex_ != -1)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneIndex()"
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
            "void Foam::solidGeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneIndex()"
            "const"
        )   << "Face zone name " << zoneName
            << " not found.  Please check your zone definition." << nl
            << "Current faceZones are:" << mesh.faceZones().names()
            << abort(FatalError);
    }

    zoneIndex_ = zone.index();
}


void Foam::solidGeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZone() const
{
    if (globalPolyPatchFaceZonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZone() const"
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
		        
}


Foam::label Foam::solidGeneralContactFvPatchVectorField::globalPolyPatchFaceZoneIndex() const
{
    if (zoneIndex_ == -1)
    {
        calcGlobalPolyPatchFaceZoneIndex();
    }

    return zoneIndex_;
}


const Foam::globalPolyPatch&
Foam::solidGeneralContactFvPatchVectorField::globalPolyPatchFaceZone() const
{
    if (!globalPolyPatchFaceZonePtr_)
    {		
		calcGlobalPolyPatchFaceZone();
    }
		
	return *globalPolyPatchFaceZonePtr_;
}


Foam::globalPolyPatch& Foam::solidGeneralContactFvPatchVectorField::globalPolyPatchFaceZone()
{

    if (!globalPolyPatchFaceZonePtr_)
    {		
		calcGlobalPolyPatchFaceZone();
    }
	
	return *globalPolyPatchFaceZonePtr_;
}


const Foam::globalPolyPatch&
Foam::solidGeneralContactFvPatchVectorField::targetGlobalPolyPatchFaceZone
(
    const label targetI
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
            field.boundaryField()[targetPatchIndices()[targetI]]
        );

    return shadowPatchField.globalPolyPatchFaceZone();
}


Foam::globalPolyPatch&
Foam::solidGeneralContactFvPatchVectorField::targetGlobalPolyPatchFaceZone
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
    solidGeneralContactFvPatchVectorField& shadowPatchField =
        const_cast<solidGeneralContactFvPatchVectorField&>
        (
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[targetPatchIndices()[targetI]]
            )
        );

    return shadowPatchField.globalPolyPatchFaceZone();
}


void Foam::solidGeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneToZones() const
{
    // Create zone-to-zone interpolation
    if (!zoneToZones_.empty())
    {
        FatalErrorIn
        (
            "void solidGeneralContactFvPatchVectorField::calcGlobalPolyPatchFaceZoneToZones()"
            "const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    zoneToZones_.setSize(targetPatchNames().size());

    const boolList& locMaster = localMaster();

    forAll(zoneToZones_, targetI)
    {
        // Only the local slave creates the interpolator
        if (locMaster[targetI])
        {
            zoneToZones_.set
                (
                    targetI,
                    new newGgiStandAlonePatchInterpolation
                    (
                        //targetGlobalPolyPatchFaceZone(targetI).globalPatch(), // master
                        globalPolyPatchFaceZone().globalPatch(), // master
						targetGlobalPolyPatchFaceZone(targetI).globalPatch(), // slave
                        tensorField(0),
                        tensorField(0),
                        vectorField(0), // Slave-to-master separation.
                        true,           // global data
                        0,              // Non-overlapping face tolerances
                        0,              //
						// Do not rescale weighting factors, as it is wrong on
						// partially covered faces						
						false,
						quickReject_,
						regionOfInterest_
                    )
                );
				
				
		
		// Check which point distance calculation method to use
            const Switch useNewPointDistanceMethod =
                dict_.lookupOrDefault<Switch>
                (
                    "useNewPointDistanceMethod", false
                );

            Info<< "    " << type() << ": " << patch().name() << nl
                << "        useNewPointDistanceMethod: "
                << useNewPointDistanceMethod
                << endl;

            zoneToZones_[targetI].useNewPointDistanceMethod() =
                useNewPointDistanceMethod;

            // Check if the projectPointsToPatchBoundary switch is set
            const Switch projectPointsToPatchBoundary =
                dict_.lookupOrDefault<Switch>
                (
                    "projectPointsToPatchBoundary",
                    false
                );

            Info<< "        projectPointsToPatchBoundary: "
                << projectPointsToPatchBoundary
                << endl;

            zoneToZones_[targetI].projectPointsToPatchBoundary() =
                projectPointsToPatchBoundary;

            if (dict_.found("checkPointDistanceOrientations"))
            {
                const Switch checkPointDistanceOrientations =
                    Switch(dict_.lookup("checkPointDistanceOrientations"));

                Info<< "        checkPointDistanceOrientations: "
                    << checkPointDistanceOrientations
                    << endl;

                zoneToZones_[targetI].checkPointDistanceOrientations() =
                    checkPointDistanceOrientations;
            }

            // Check if the usePrevCandidateMasterNeighbors switch is set
            const Switch usePrevCandidateMasterNeighbors =
                dict_.lookupOrDefault<Switch>
                (
                    "usePrevCandidateMasterNeighbors",
                    false
                );

            Info<< "        usePrevCandidateMasterNeighbors: "
                << usePrevCandidateMasterNeighbors
                << endl;

            zoneToZones_[targetI].usePrevCandidateMasterNeighbors() =
                usePrevCandidateMasterNeighbors;
		
		
        }
    }
}


const Foam::newGgiStandAlonePatchInterpolation&
Foam::solidGeneralContactFvPatchVectorField::zoneToZone
(
    const label targetI
) const
{
	
    if (!localMaster()[targetI])
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
Foam::solidGeneralContactFvPatchVectorField::zoneToZone(const label targetI)
{
    if (!localMaster()[targetI])
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


Foam::label Foam::solidGeneralContactFvPatchVectorField::findShadowID
(
    const label patchID
) const
{
    label targetI = -1;

    const labelList shadowIDs = targetPatchIndices();

    forAll(shadowIDs, I)
    {
        if (patchID == shadowIDs[I])
        {
            targetI = I;
            break;
        }
    }

    if (targetI == -1)
    {
        FatalErrorIn("findShadowID(const label patchID)")
            << "shadow patch not found!" << abort(FatalError);
    }

    return targetI;
}


void Foam::solidGeneralContactFvPatchVectorField::makeCurPatchTractions() const
{
    if (curPatchTractionPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
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

    const boolList& locMaster = localMaster();

    forAll(locMaster, targetI)
    {
        vectorField curPatchSlip(patch().size(), vector::zero);

        // Calculate slip
        if (locMaster[targetI])
        {
            curPatchSlip = frictionModel(targetI).slip();
        }
        else
        {
            const solidGeneralContactFvPatchVectorField& shadowPatchField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[targetPatchIndices()[targetI]]
            );

            const label locShadowID =
                shadowPatchField.findShadowID(patch().index());

            vectorField shadowPatchSlip =
                shadowPatchField.frictionModel(locShadowID).slip();

            vectorField shadowZoneSlip = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(shadowPatchSlip);
			

            // Interpolate from shadow to the current patch
            // Face-to-face
			

            vectorField curZoneSlip =
                shadowPatchField.zoneToZone(locShadowID).slaveToMaster
                (
                    shadowZoneSlip
                );

            curPatchSlip = globalPolyPatchFaceZone().globalFaceToPatch(curZoneSlip);
			
        }

        // Heat flux rate: rate of dissipated frictional energy
        // The dot product of the traction vectors and the slip vectors
        // gives the dissipated frictional energy rate per unit area, which
        // is always positive
        Qc += mag(curTraction & (curPatchSlip/deltaT));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localMasterPtr_(NULL),
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
	quickReject_(Foam::newGgiInterpolation::AABB),
	regionOfInterestTopCorner_(vector::max),
    regionOfInterestBottomCorner_(vector::min),
    regionOfInterest_(vector::min, vector::max),
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
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        ) << endl;
    }
}


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localMasterPtr_(NULL),
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
	quickReject_(ptf.quickReject_),
	regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
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

    if (ptf.firstPatchInListPtr_)
    {
        firstPatchInListPtr_ = new bool(*ptf.firstPatchInListPtr_);
    }

    if (ptf.firstPatchInListIndexPtr_)
    {
        firstPatchInListIndexPtr_ = new label(*ptf.firstPatchInListIndexPtr_);
    }

    if (ptf.localMasterPtr_)
    {
        localMasterPtr_ = new boolList(*ptf.localMasterPtr_);
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
	

    if (!ptf.zoneToZones_.empty())
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


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localMasterPtr_(NULL),
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
	quickReject_
    (
        newGgiInterpolation::quickRejectNames_
        [
            dict.lookupOrDefault<word>("quickReject", "AABB")
        ]
    ),
	regionOfInterestTopCorner_
    (
        dict.lookupOrDefault<vector>
        (
            "regionOfInterestTopCorner",
            vector::max
        )
    ),
    regionOfInterestBottomCorner_
    (
        dict.lookupOrDefault<vector>
        (
            "regionOfInterestBottomCorner",
            vector::min
        )
    ),
    regionOfInterest_
    (
        boundBox
        (
            regionOfInterestBottomCorner_,
            regionOfInterestTopCorner_
        )
//        dict.lookupOrDefault<boundBox>
//        (
//            "regionOfInterest",
//            boundBox(vector::min, vector::max)
//        )
    ),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(0.0)
{
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


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf
)
:
    solidTractionFvPatchVectorField(ptf),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localMasterPtr_(NULL),
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
	quickReject_(ptf.quickReject_),
	regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
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
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf"
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

    if (ptf.localMasterPtr_)
    {
        localMasterPtr_ = new boolList(*ptf.localMasterPtr_);
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


    if (!ptf.zoneToZones_.empty())
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


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
    firstPatchInListPtr_(NULL),
    firstPatchInListIndexPtr_(NULL),
    localMasterPtr_(NULL),
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
	quickReject_(ptf.quickReject_),
	regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
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
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
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

    if (ptf.localMasterPtr_)
    {
        localMasterPtr_ = new boolList(*ptf.localMasterPtr_);
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


    if (!ptf.zoneToZones_.empty())
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


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


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

    deleteDemandDrivenData(firstPatchInListPtr_);
    deleteDemandDrivenData(firstPatchInListIndexPtr_);
    deleteDemandDrivenData(localMasterPtr_);
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


bool Foam::solidGeneralContactFvPatchVectorField::firstPatchInList() const
{
    if (!firstPatchInListPtr_)
    {
        calcFirstPatchInList();
    }

    return *firstPatchInListPtr_;
}


Foam::label Foam::solidGeneralContactFvPatchVectorField::firstPatchInListIndex()
const
{
    if (!firstPatchInListIndexPtr_)
    {
        calcFirstPatchInListIndex();
    }

    return *firstPatchInListIndexPtr_;
}


const Foam::List<Foam::word>&
Foam::solidGeneralContactFvPatchVectorField::targetPatchNames() const
{
    if (!targetPatchNamesPtr_)
    {
        calcTargetPatchNames();
    }

    return *targetPatchNamesPtr_;
}


const Foam::List<Foam::label>&
Foam::solidGeneralContactFvPatchVectorField::targetPatchIndices() const
{
    if (!targetPatchIndicesPtr_)
    {
        calcTargetPatchNames();
    }

    return *targetPatchIndicesPtr_;
}


const Foam::List<Foam::word>&
Foam::solidGeneralContactFvPatchVectorField::targetGlobalPolyPatchFaceZoneNames() const
{
    if (!targetGlobalPolyPatchFaceZoneNamesPtr_)
    {
        calcTargetGlobalPolyPatchFaceZoneNames();
    }

    return *targetGlobalPolyPatchFaceZoneNamesPtr_;
}


const Foam::List<Foam::label>&
Foam::solidGeneralContactFvPatchVectorField::targetGlobalPolyPatchFaceZoneIndices() const
{
    if (!targetGlobalPolyPatchFaceZoneIndicesPtr_)
    {
        calcTargetGlobalPolyPatchFaceZoneNames();
    }

    return *targetGlobalPolyPatchFaceZoneIndicesPtr_;
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
}


void Foam::solidGeneralContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    boolList activeContactPairs(targetPatchNames().size(), true);  //true

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
				
				// Force N^2 contact search at least once per time-step
                zoneToZone(activePatchI).clearPrevCandidateMasterNeighbors();
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
            if (localMaster()[slaveI])
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
        const boolList& locMaster = localMaster();

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
			
            // Perform quick check to find potential contacting pairs
            // The quick check is based on the bounding box (BB) of the contact
            // pairs: if the BBs of pair intersect then we will designate the
            // pair as active.
			
			
			
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
			

            // Call normal and frction contact models for active contacting
            // pairs
            // Accumulate contact force contributions for all active contact
            // pairs

            if (activeContactPairs[activePatchI])
            {
                if (locMaster[activePatchI])
                {
					// Correct normal and friction contact models for the
                    // current contact pair

                    // Calculate the slave patch face unit normals as they are
                    // units by both the normal and friction models
                    const vectorField targetPatchFaceNormals = 
										targetGlobalPolyPatchFaceZone(activePatchI).globalFaceToPatch
										(
										targetGlobalPolyPatchFaceZone(activePatchI).globalPatch().faceNormals()
										);
					
					
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
					
					
					
                    // Master zone DD
                    const vectorField zoneDD = globalPolyPatchFaceZone().patchFaceToGlobal(patchDD);
					
					
                    // Master patch DD interpolated to the slave patch
                    const vectorField patchDDInterpToShadowPatch =
						targetGlobalPolyPatchFaceZone(activePatchI).globalFaceToPatch
										(
										zoneToZone(activePatchI).masterToSlave(zoneDD)()
										);
					
					 normalModel(activePatchI).correct
                     (
                         targetPatchFaceNormals,
                         //zoneToZone(activePatchI),
						 targetGlobalPolyPatchFaceZone(activePatchI).globalPointToPatch
						 (
						 zoneToZone(activePatchI).slavePointDistanceToIntersection()
						 ),
                         targetPatchDD,
                         patchDDInterpToShadowPatch
                     );
					 
					 
                    frictionModel(activePatchI).correct
                    (
                        normalModel(activePatchI).slavePressure(),
                        targetPatchFaceNormals,
                        normalModel(activePatchI).areaInContact(),
                        targetPatchDD,
                        patchDDInterpToShadowPatch
                    );
					
					
				
				vectorField targetPatchTraction =
                        - frictionModel(activePatchI).slaveTractionForMaster()
                        - normalModel(activePatchI).slavePressure();

                    vectorField shadowZoneTraction = targetGlobalPolyPatchFaceZone(activePatchI).patchFaceToGlobal(targetPatchTraction);
                        
					
                    // Face-to-face
                    vectorField masterZoneTraction = 
					globalPolyPatchFaceZone().globalFaceToPatch
					(
					zoneToZone(activePatchI).slaveToMaster(shadowZoneTraction)()
					);
					
                    // We store master patch traction as thermalGeneralContact
                    // uses it
                    curPatchTractions(activePatchI) =
						globalPolyPatchFaceZone().globalFaceToPatch(masterZoneTraction);

                    curPatchTraction += curPatchTractions(activePatchI);
					
				}
                else // local slave
                {
					
					// Get traction from local slave
	
                    const volVectorField& field =
                        db().lookupObject<volVectorField>
                        (
                            dimensionedInternalField().name()
                        );
						
					// Lookup the master patch corresponding to the current slave patch
                    const solidGeneralContactFvPatchVectorField&
                        localMasterField =
                        refCast<const solidGeneralContactFvPatchVectorField>
                        (
                            field.boundaryField()
                            [
                                targetPatchIndices()[activePatchI]
                            ]
                        );

                    const label slaveShadowI =
                        localMasterField.findShadowID(patch().index());
						
					// Accumulate traction

                    curPatchTractions(activePatchI) =
                        localMasterField.frictionModel(slaveShadowI).slaveTraction()
                        + localMasterField.normalModel(slaveShadowI).slavePressure();
						
					curPatchTraction += curPatchTractions(activePatchI);
						
					
					
                }
            } // if contact pair is active
        } // forAll contact pairs

        // Set master gradient based on accumulated traction
        traction() = curPatchTraction;
    }

    solidTractionFvPatchVectorField::updateCoeffs();
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
    const label targetI
) const
{
    if (!QcsPtr_)
    {
        calcQcs();
    }

    return (*QcsPtr_)[targetI];
}


void Foam::solidGeneralContactFvPatchVectorField::calcQcs() const
{
    const boolList& locMaster = localMaster();

    QcsPtr_ = new List<scalarField>(locMaster.size());

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

    forAll(locMaster, targetI)
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
        if (locMaster[targetI])
        {
            curPatchSlip = frictionModel(targetI).slip();
        }
        else
        {
            const solidGeneralContactFvPatchVectorField& shadowPatchField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[targetPatchIndices()[targetI]]
            );

            const label locShadowID =
            shadowPatchField.findShadowID(patch().index());

            vectorField shadowPatchSlip =
            shadowPatchField.frictionModel(locShadowID).slip();

            vectorField shadowZoneSlip = targetGlobalPolyPatchFaceZone(targetI).patchFaceToGlobal(shadowPatchSlip);


            // Interpolate from shadow to the current patch
            // Face-to-face
			
			
            vectorField curZoneSlip =
            shadowPatchField.zoneToZone(locShadowID).slaveToMaster
            (
                shadowZoneSlip
            );

            curPatchSlip = globalPolyPatchFaceZone().globalFaceToPatch(curZoneSlip);

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

    //bbOffset_ = 5.0*returnReduce(minDim, minOp<scalar>()); //oRG
	bbOffset_ = 1.25*returnReduce(minDim, minOp<scalar>());

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
Foam::solidGeneralContactFvPatchVectorField::curPatchTractions
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


// Write
void Foam::solidGeneralContactFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("rigidMaster")
        << rigidMaster_ << token::END_STATEMENT << nl;

    // Write the dict from the first contact model

    const label targetI = 0;

    if (localMaster()[targetI])
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

        const solidGeneralContactFvPatchVectorField& localMasterField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()
                [
                    targetPatchIndices()[targetI]
                ]
            );

        const label localMasterID =
            localMasterField.findShadowID(patch().index());

        os.writeKeyword("normalContactModel")
            << localMasterField.normalModel(localMasterID).type()
            << token::END_STATEMENT << nl;
        localMasterField.normalModel(localMasterID).writeDict(os);

        os.writeKeyword("frictionContactModel")
            << localMasterField.frictionModel(localMasterID).type()
            << token::END_STATEMENT << nl;
        localMasterField.frictionModel(localMasterID).writeDict(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField, solidGeneralContactFvPatchVectorField
    )
}


// ************************************************************************* //
