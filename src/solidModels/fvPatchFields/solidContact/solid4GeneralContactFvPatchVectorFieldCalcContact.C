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

InClass
    solid4GeneralContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solid4GeneralContactFvPatchVectorField.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::solid4GeneralContactFvPatchVectorField::moveZonesToDeformedConfiguration()
{
	Info<<"IN -- moveZonesToDeformedConfiguration() line:"<<__LINE__<<endl;
    // Only the master moves the zones
    if (!firstPatchInList())
    {
        return;
    }

    // Method
    // We will interpolate the patch face displacements to the patch vertices
    // and then add these vertex/point displacements to the initial patch
    // points
    // We need to take care in parallel, and also realise that the solidModel
    // might have a moving or stationary mesh

    // Assemble the zone face displacement field to move the zones

    // First we will move the master zone

    vectorField zoneD(zone().globalPatch().size(), vector::zero);

    // For a non-moving mesh, we will move the zones by the total
    // displacement, whereas for a moving mesh (updated Lagrangian), we will
    // move the zones by the displacement increment

    if (movingMesh())
    {
        // Updated Lagrangian, so we will move the zones by the displacement
        // increment

        // Lookup the current total displacement field 
        const volVectorField& DD = db().lookupObject<volVectorField>("DU");
	//	const volVectorField& DD = db().lookupObject<volVectorField>("DD");

        // Take a reference to the patch face displacement increment field
        const vectorField& patchDD =
            DD.boundaryField()[patch().index()];

        zoneD = zone().patchFaceToGlobal(patchDD);
    }
    else
    {
        // Non-moving mesh: we will move the zones by the total displacement

        // Lookup the current total displacement field
        const volVectorField& D = db().lookupObject<volVectorField>("U");
	//	const volVectorField& D = db().lookupObject<volVectorField>("D");

        // Take a reference to the patch face total displacement field
        const vectorField& patchD =
            D.boundaryField()[patch().index()];

        Info<<"patchD IN -- moveZonesToDeformedConfiguration(): "<<patchD<<endl;
		
		zoneD = zone().patchFaceToGlobal(patchD);
    }

    // Interpolate the zone face field to the zone points
    const pointField zonePointD =
        zone().interpolator().faceToPointInterpolate(zoneD);

    // The zone deformed points are the initial position plus the
    // displacement
    const pointField zoneNewPoints =
        zone().patchPointToGlobal
        (
            patch().patch().localPoints()
        )
      + zonePointD;

    // Remove zone weights
    zone().movePoints(zoneNewPoints);

    // We need to use const_cast to move the standAlonePatch points as the
    // movePoints function only clears weights
    // Also, be careful to move the points as opposed to the localPoints
    const_cast<pointField&>(zone().globalPatch().points()) = zoneNewPoints;

    // Secondly we will move the slave zones

    forAll(slavePatchNames(), shadPatchI)
    {
        vectorField slaveZoneD
        (
            slaveZones()[shadPatchI].globalPatch().size(),
            vector::zero
        );

        // For a non-moving mesh, we will move the zones by the total
        // displacement, whereas for a moving mesh (updated Lagrangian), we will
        // move the zones by the displacement increment

        if (movingMesh())
        {
            // Updated Lagrangian, so we will move the zones by the displacement
            // increment

            // Lookup the current total displacement field
            const volVectorField& DD = db().lookupObject<volVectorField>("DU");
		//	const volVectorField& DD = db().lookupObject<volVectorField>("DD");

            // Take a reference to the patch face displacement increment field
            const vectorField& slavePatchDD =
                DD.boundaryField()[slavePatchIndices()[shadPatchI]];

            slaveZoneD = slaveZones()[shadPatchI].patchFaceToGlobal
            (
                slavePatchDD
            );
        }
        else
        {
            // Non-moving mesh: we will move the zones by the total displacement

            // Lookup the current total displacement field
            const volVectorField& D = db().lookupObject<volVectorField>("U");
		//	const volVectorField& D = db().lookupObject<volVectorField>("D");
			Info<<"D.size() internal field size in moveZonesToDeformedConfiguration(): "<<D.size()<<endl;
			
            // Take a reference to the patch face total displacement field
            const vectorField& slavePatchD =
                D.boundaryField()[slavePatchIndices()[shadPatchI]];

            slaveZoneD = slaveZones()[shadPatchI].patchFaceToGlobal
            (
                slavePatchD
            );
        }

        // Interpolate the zone face field to the zone points
        const pointField slaveZonePointD =
            slaveZones()[shadPatchI].interpolator().faceToPointInterpolate
            (
                slaveZoneD
            );

        // The zone deformed points are the initial position plus the
        // displacement
        const pointField slaveZoneNewPoints =
            slaveZones()[shadPatchI].patchPointToGlobal
            (
                slaveZones()[shadPatchI].patch().localPoints()
            )
          + slaveZonePointD;

        // Remove zone weights
        slaveZones()[shadPatchI].movePoints(slaveZoneNewPoints);

        // We need to use const_cast to move the standAlonePatch points as the
        // movePoints function only clears weights
        // Also, be careful to move the points are opposed to the localPoints
        const_cast<pointField&>
        (
            slaveZones()[shadPatchI].globalPatch().points()
        ) = slaveZoneNewPoints;
    }
	Info<<"IN -- moveZonesToDeformedConfiguration() line:"<<__LINE__<<endl;
}


void Foam::solid4GeneralContactFvPatchVectorField::calcZone() const
{
    if (debug)
    {
        InfoIn("void Foam::solid4GeneralContactFvPatchVectorField::calcZone() const")
            << patch().name() << " : making the zone" << endl;
    }
	
	// Remove this check for now
	
    if (!currentMaster())
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::calcZone() const"
        )   << "Trying to create zone on a slave" << abort(FatalError);
    }
	

    if (zonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::calcZone() const"
        )   << "pointer already set" << abort(FatalError);
    }

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    zonePtr_ = new globalPolyPatch
    (
        patch().name(),
        patch().boundaryMesh().mesh()
    );
}


void Foam::solid4GeneralContactFvPatchVectorField::calcSlaveZones() const
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::calcSlaveZones() const"
        )   << patch().name() << " : making the slave zones" << endl;
    }

    if (!currentMaster())
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::"
            "calcSlaveZones() const"
        )   << "Trying to create slave zones on a slave" << abort(FatalError);
    }

    if (!slaveZones_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::"
            "calcSlaveZones() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const wordList& shadPatchNames = slavePatchNames();

    slaveZones_.setSize(shadPatchNames.size());

    forAll(slaveZones_, shadPatchI)
    {
        // Note: the main mesh will either be in the initial configuration or
        // the updated configuration
        slaveZones_.set
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


void Foam::solid4GeneralContactFvPatchVectorField::calcZoneToZones() const
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::calcZoneToZones() const"
        )   << patch().name() << " : making the zoneToZone" << endl;
    }

    // Create zone-to-zone interpolation
    if (!zoneToZones_.empty())
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::calcZoneToZones() const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }
	
	Info<< "patch().name() in calcZoneToZones() "<<patch().name()<<endl;
	Info<< "patch().index() in calcZoneToZones() "<<patch().index()<<endl;
	Info<<"In calcZoneToZones() line: "<<__LINE__<<endl;
    
	// Check master and slave patch
    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );
	
	
	const boolList& locSlave = localSlave();

    zoneToZones_.setSize(slavePatchNames().size());

    forAll (zoneToZones_, shadPatchI)
    {
        
		const solid4GeneralContactFvPatchVectorField& slavePatchField =
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()[slavePatchIndices()[shadPatchI]]
            );

        if (currentMaster())
        {
            if (slavePatchField.master() == true)
            {
                FatalErrorIn
                (
                    "void solidContactFvPatchScalarField::"
                    "calcZoneToZones() const"
                )   << "There are two master patches!" << abort(FatalError);
            }
        }
        else
        {
            if (slavePatchField.master() == false)
            {
                FatalErrorIn
                (
                    "void solidContactFvPatchScalarField::"
                    "calcZoneToZones() const"
                )   << "There is no master patch!" << abort(FatalError);
            }
        }
		
        
		if (locSlave[shadPatchI]) //if (firstPatchInList())
        {
			Info<<"In calcZoneToZones() line: "<<__LINE__<<endl;
            // Create interpolation for patches
            zoneToZones_.set
            (
                shadPatchI,
                new newGgiStandAlonePatchInterpolation
                (
                    zone().globalPatch(),
                    slaveZones()[shadPatchI].globalPatch(),
                    tensorField(0),
                    tensorField(0),
                    vectorField(0), // Slave-to-master separation. Bug fix
                    true,           // global data
                    0,              // Master non-overlapping face tolerances
                    0,              // Slave non-overlapping face tolerances
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

            zoneToZones_[shadPatchI].useNewPointDistanceMethod() =
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

            zoneToZones_[shadPatchI].projectPointsToPatchBoundary() =
                projectPointsToPatchBoundary;

            if (dict_.found("checkPointDistanceOrientations"))
            {
                const Switch checkPointDistanceOrientations =
                    Switch(dict_.lookup("checkPointDistanceOrientations"));

                Info<< "        checkPointDistanceOrientations: "
                    << checkPointDistanceOrientations
                    << endl;

                zoneToZones_[shadPatchI].checkPointDistanceOrientations() =
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

            zoneToZones_[shadPatchI].usePrevCandidateMasterNeighbors() =
                usePrevCandidateMasterNeighbors;
        }
        else
        {
            FatalErrorIn
            (
                "void solid4GeneralContactFvPatchVectorField::"
                "calcZoneToZones() const"
            )   << "Attempting to create GGIInterpolation on a slave patch"
                << abort(FatalError);
        }
    }
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

    contactPerSlave_.setSize(slavePatchNames().size());

    forAll(contactPerSlave_, i)
    {
        contactPerSlave_.set
        (
            i,
            new scalarField(patch().size(), 0.0)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::globalPolyPatch&
Foam::solid4GeneralContactFvPatchVectorField::zone() const
{
	Info<< "patch().name() in zone() "<<patch().name()<<endl;
	Info<< "patch().index() in zone() "<<patch().index()<<endl;
	Info<<"In zone() line:"<<__LINE__<<endl;
    
	if (currentMaster())
    {
		Info<<"In zone() line:"<<__LINE__<<endl;
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }
    
	else
    {
		Info<<"In zone() line:"<<__LINE__<<endl;
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        const solid4GeneralContactFvPatchVectorField& slavePatchField =
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()[slavePatchIndices()[0]]
            );

        return slavePatchField.zone();
    }
	
}


Foam::globalPolyPatch& Foam::solid4GeneralContactFvPatchVectorField::zone()
{
	Info<< "patch().name() in zone() "<<patch().name()<<endl;
	Info<< "patch().index() in zone() "<<patch().index()<<endl;
	Info<<"In zone() line:"<<__LINE__<<endl;
			
    if (currentMaster())
    {
		Info<<"In zone() line:"<<__LINE__<<endl;
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    
	}	
    else
    {
		Info<<"In zone() line:"<<__LINE__<<endl;
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        solid4GeneralContactFvPatchVectorField& slavePatchField =
            const_cast<solid4GeneralContactFvPatchVectorField&>
            (
                refCast<const solid4GeneralContactFvPatchVectorField>
                (
                    field.boundaryField()[slavePatchIndices()[0]]
                )
            );
		
        return slavePatchField.zone();
    }	
}


const Foam::PtrList<Foam::globalPolyPatch>&
Foam::solid4GeneralContactFvPatchVectorField::slaveZones() const
{
    if (currentMaster())
    {
        if (slaveZones_.empty())
        {
            calcSlaveZones();
        }

        return slaveZones_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        const solid4GeneralContactFvPatchVectorField& slavePatchField =
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()[slavePatchIndices()[0]]
            );

        return slavePatchField.slaveZones();
    }
}


Foam::PtrList<Foam::globalPolyPatch>&
Foam::solid4GeneralContactFvPatchVectorField::slaveZones()
{
    if (currentMaster())
    {
        if (slaveZones_.empty())
        {
            calcSlaveZones();
        }
		
        return slaveZones_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        solid4GeneralContactFvPatchVectorField& slavePatchField =
            const_cast<solid4GeneralContactFvPatchVectorField&>
            (
                refCast<const solid4GeneralContactFvPatchVectorField>
                (
                    field.boundaryField()[slavePatchIndices()[0]]
                )
            );

        return slavePatchField.slaveZones();
    }
}


const Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
Foam::solid4GeneralContactFvPatchVectorField::zoneToZones() const
{
	Info<<"IN -- zoneToZones() line:"<<__LINE__<<endl;
    if (currentMaster())
    {
        if (zoneToZones_.empty())
        {
            calcZoneToZones();
        }

        return zoneToZones_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        const solid4GeneralContactFvPatchVectorField& slavePatchField =
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()[slavePatchIndices()[0]]
            );

        return slavePatchField.zoneToZones();
    }
}


Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
Foam::solid4GeneralContactFvPatchVectorField::zoneToZones()
{
    Info<<"IN -- zoneToZones() line:"<<__LINE__<<endl;
    if (currentMaster())
    {
        if (zoneToZones_.empty())
        {
            calcZoneToZones();
        }

        return zoneToZones_;
    }
    else
    {
        // We will const_cast the slave patch so we can delete the weights when
        // the zones move
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        solid4GeneralContactFvPatchVectorField& slavePatchField =
            const_cast<solid4GeneralContactFvPatchVectorField&>
            (
                refCast<const solid4GeneralContactFvPatchVectorField>
                (
                    field.boundaryField()[slavePatchIndices()[0]]
                )
            );

        return slavePatchField.zoneToZones();
    }

}


const Foam::newGgiStandAlonePatchInterpolation&
Foam::solid4GeneralContactFvPatchVectorField::zoneToZoneForThisSlave() const
{
    if (currentMaster())
    {
        FatalErrorIn
        (
            "const Foam::newGgiStandAlonePatchInterpolation&"
            " Foam::solid4GeneralContactFvPatchVectorField::"
            "zoneToZoneForThisSlave() const"
        )   << "The master patch is not allowed to call this function"
            << abort(FatalError);
    }

    // The master may have multiple slaves so we need to find which zoneToZone
    // corresponds to the current slave patch
    const wordList& shadPatchNames = slavePatchField().slavePatchNames();
    label masterSlaveID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterSlaveID = shadPatchI;
            break;
        }
    }

    if (masterSlaveID == -1)
    {
        FatalErrorIn
        (
            "const Foam::newGgiStandAlonePatchInterpolation&"
            " Foam::solid4GeneralContactFvPatchVectorField::"
            "zoneToZoneForThisSlave() const"
        )   << "Something went wrong when looking for the slavePatch"
            << abort(FatalError);
    }

    // Return the zoneToZone between the master and the current patch
    return zoneToZones()[masterSlaveID];
}


const Foam::globalPolyPatch&
Foam::solid4GeneralContactFvPatchVectorField::zoneForThisSlave() const
{
    if (currentMaster())
    {
        FatalErrorIn
        (
            "const Foam::globalPolyPatch&"
            " Foam::solid4GeneralContactFvPatchVectorField::zoneForThisSlave() const"
        )   << "The master patch is not allowed to call this function"
            << abort(FatalError);
    }

    // The master may have multiple slaves so we need to find which slaveZone
    // corresponds to the current slave patch
    const wordList& shadPatchNames = slavePatchField().slavePatchNames();
    label masterSlaveID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterSlaveID = shadPatchI;
            break;
        }
    }

    if (masterSlaveID == -1)
    {
        FatalErrorIn
        (
            "const Foam::globalPolyPatch&"
            " Foam::solid4GeneralContactFvPatchVectorField::zoneForThisSlave() const"
        )   << "Something went wrong when looking for the slavePatch"
            << abort(FatalError);
    }

    // Return the zoneToZone between the master and the current patch
    return slaveZones()[masterSlaveID];
}


//*************** based on solidGeneral*****************
/*
// Private GeneralMember functions
void Foam::solid4GeneralContactFvPatchVectorField::calcslaveGPatchNames() const
{
    if (slaveGPatchNamesPtr_ || slaveGPatchIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solid4GeneralContactFvPatchVectorField::"
            "calcslaveGPatchNames() const"
        )   << "slaveGPatchNames_ or slaveGPatchIndices_ already set"
            << abort(FatalError);
    }

    // Add each solid4GeneralContact patch in the order of increasing patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    // Count slaveG patches

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
	
    slaveGPatchNamesPtr_ = new wordList(nShadPatches);
    wordList& slaveGPatchNames = *slaveGPatchNamesPtr_;
	Info<<"In calcslaveGPatchNames():"<<__LINE__<<endl;
//	Info<<"What is *slaveGPatchNamesPtr_? in calcslaveGPatchNames():"<<*slaveGPatchNamesPtr_<<endl;

    slaveGPatchIndicesPtr_ = new labelList(nShadPatches);
    labelList& slaveGPatchIndices = *slaveGPatchIndicesPtr_;
//	Info<<"In calcslaveGPatchNames():"<<__LINE__<<endl;
	
    // Record slaveG patch names

    label slaveGI = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solid4GeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
            slaveGPatchNames[slaveGI] = patch().boundaryMesh()[patchI].name();

            slaveGPatchIndices[slaveGI++] = patchI;
        }
    }
}

void Foam::solid4GeneralContactFvPatchVectorField::calcslaveGZoneNames() const
{
    if (slaveGZoneNamesPtr_ || slaveGZoneIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solid4GeneralContactFvPatchVectorField::"
            "calcslaveGZoneNames() const"
        )   << "slaveGZoneNames_ or slaveGZoneIndices_ already set"
            << abort(FatalError);
    }

    const wordList& shadNames = slaveGPatchNames();

    slaveGZoneNamesPtr_ = new wordList(shadNames.size());
    wordList& slaveGZoneNames = *slaveGZoneNamesPtr_;

    slaveGZoneIndicesPtr_ = new labelList(shadNames.size());
    labelList& slaveGZoneIndices = *slaveGZoneIndicesPtr_;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    forAll(shadNames, slaveGI)
    {
        word zoneName = shadNames[slaveGI] + "FaceZone";

        faceZoneID zone(zoneName, mesh.faceZones());

        if (!zone.active())
        {
            FatalErrorIn("solid4GeneralContactFvPatchVectorField")
                << "Face zone name " << zoneName
                << " not found.  Please check your zone definition."
                << abort(FatalError);
        }

        slaveGZoneNames[slaveGI] = zoneName;

        slaveGZoneIndices[slaveGI] = zone.index();
    }
}
		
//Public:
//Access GeneralMember

//Step1 - 
const Foam::List<Foam::word>&
Foam::solid4GeneralContactFvPatchVectorField::slaveGPatchNames() const
{
    if (!slaveGPatchNamesPtr_)
    {
        calcslaveGPatchNames();
    }

    return *slaveGPatchNamesPtr_;
}

//Step3 - 
const Foam::List<Foam::label>&
Foam::solid4GeneralContactFvPatchVectorField::slaveGPatchIndices() const
{
	Info<<"In slaveGPatchIndices() line "<<__LINE__<<endl;
    if (!slaveGPatchIndicesPtr_)
    {
		Info<<"In slaveGPatchIndices() line "<<__LINE__<<endl;
        calcslaveGPatchNames();
    }
Info<<"*slaveGPatchIndicesPtr_ in slaveGPatchIndices() "<<*slaveGPatchIndicesPtr_<<endl;
    return *slaveGPatchIndicesPtr_;
}
        
const Foam::List<Foam::word>&
Foam::solid4GeneralContactFvPatchVectorField::slaveGZoneNames() const
{
    if (!slaveGZoneNamesPtr_)
    {
        calcslaveGZoneNames();
    }

    return *slaveGZoneNamesPtr_;
}

const Foam::List<Foam::label>&
Foam::solid4GeneralContactFvPatchVectorField::slaveGZoneIndices() const
{
    if (!slaveGZoneIndicesPtr_)
    {
        calcslaveGZoneNames();
    }

    return *slaveGZoneIndicesPtr_;
}

			
//Step2's - 
const Foam::standAlonePatch&
Foam::solid4GeneralContactFvPatchVectorField::slaveGZone
(
    const label slaveGI
) const
{
    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    const solid4GeneralContactFvPatchVectorField& slaveGPatchField =
        refCast<const solid4GeneralContactFvPatchVectorField>
        (
            field.boundaryField()[slaveGPatchIndices()[slaveGI]]
        );

    return slaveGPatchField.zone();
}


Foam::standAlonePatch&
Foam::solid4GeneralContactFvPatchVectorField::slaveGZone
(
    const label slaveGI
)
{
    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    // Const cast away the const-ness
    solid4GeneralContactFvPatchVectorField& slaveGPatchField =
        const_cast<solid4GeneralContactFvPatchVectorField&>
        (
            refCast<const solid4GeneralContactFvPatchVectorField>
            (
                field.boundaryField()[slaveGPatchIndices()[slaveGI]]
            )
        );

    return slaveGPatchField.zone();
}

*/
Foam::label Foam::solid4GeneralContactFvPatchVectorField::findSlaveID
(
    const label patchID
) const
{
    label shadPatchI = -1;

    const labelList slaveIDs = slavePatchIndices();

    forAll(slaveIDs, I)
    {
        if (patchID == slaveIDs[I])
        {
            shadPatchI = I;
            break;
        }
    }

    if (shadPatchI == -1)
    {
        FatalErrorIn("findSlaveID(const label patchID)")
            << "slave patch not found!" << abort(FatalError);
    }

    return shadPatchI;
}

//*************** END based on solidGeneral**************


// ************************************************************************* //
