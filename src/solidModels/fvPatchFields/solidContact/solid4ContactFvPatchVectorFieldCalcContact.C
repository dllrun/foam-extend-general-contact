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
    solid4ContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solid4ContactFvPatchVectorField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::solid4ContactFvPatchVectorField::moveZonesToDeformedConfiguration()
{
    // Only the master moves the zones
    if (!master_)
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

    // Secondly we will move the shadow zones

    forAll(shadowPatchNames(), shadPatchI)
    {
        vectorField shadowZoneD
        (
            shadowZones()[shadPatchI].globalPatch().size(),
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
            const vectorField& shadowPatchDD =
                DD.boundaryField()[shadowPatchIndices()[shadPatchI]];

            shadowZoneD = shadowZones()[shadPatchI].patchFaceToGlobal
            (
                shadowPatchDD
            );
        }
        else
        {
            // Non-moving mesh: we will move the zones by the total displacement

            // Lookup the current total displacement field
            const volVectorField& D = db().lookupObject<volVectorField>("U");
		//	const volVectorField& D = db().lookupObject<volVectorField>("D");

            // Take a reference to the patch face total displacement field
            const vectorField& shadowPatchD =
                D.boundaryField()[shadowPatchIndices()[shadPatchI]];

            shadowZoneD = shadowZones()[shadPatchI].patchFaceToGlobal
            (
                shadowPatchD
            );
        }

        // Interpolate the zone face field to the zone points
        const pointField shadowZonePointD =
            shadowZones()[shadPatchI].interpolator().faceToPointInterpolate
            (
                shadowZoneD
            );

        // The zone deformed points are the initial position plus the
        // displacement
        const pointField shadowZoneNewPoints =
            shadowZones()[shadPatchI].patchPointToGlobal
            (
                shadowZones()[shadPatchI].patch().localPoints()
            )
          + shadowZonePointD;

        // Remove zone weights
        shadowZones()[shadPatchI].movePoints(shadowZoneNewPoints);

        // We need to use const_cast to move the standAlonePatch points as the
        // movePoints function only clears weights
        // Also, be careful to move the points are opposed to the localPoints
        const_cast<pointField&>
        (
            shadowZones()[shadPatchI].globalPatch().points()
        ) = shadowZoneNewPoints;
    }
}


void Foam::solid4ContactFvPatchVectorField::calcZone() const
{
    if (debug)
    {
        InfoIn("void Foam::solid4ContactFvPatchVectorField::calcZone() const")
            << patch().name() << " : making the zone" << endl;
    }

    if (!master_)
    {
        FatalErrorIn
        (
            "void Foam::solid4ContactFvPatchVectorField::calcZone() const"
        )   << "Trying to create zone on a slave" << abort(FatalError);
    }

    if (zonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solid4ContactFvPatchVectorField::calcZone() const"
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


void Foam::solid4ContactFvPatchVectorField::calcShadowZones() const
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solid4ContactFvPatchVectorField::calcShadowZones() const"
        )   << patch().name() << " : making the shadow zones" << endl;
    }

    if (!master_)
    {
        FatalErrorIn
        (
            "void Foam::solid4ContactFvPatchVectorField::"
            "calcShadowZones() const"
        )   << "Trying to create shadow zones on a slave" << abort(FatalError);
    }

    if (!shadowZones_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solid4ContactFvPatchVectorField::"
            "calcShadowZones() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const wordList& shadPatchNames = shadowPatchNames();

    shadowZones_.setSize(shadPatchNames.size());

    forAll(shadowZones_, shadPatchI)
    {
        // Note: the main mesh will either be in the initial configuration or
        // the updated configuration
        shadowZones_.set
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


void Foam::solid4ContactFvPatchVectorField::calcZoneToZones() const
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solid4ContactFvPatchVectorField::calcZoneToZones() const"
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

    // Check master and slave patch
    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    zoneToZones_.setSize(shadowPatchNames().size());

    forAll (zoneToZones_, shadPatchI)
    {
        const solid4ContactFvPatchVectorField& shadowPatchField =
            refCast<const solid4ContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadPatchI]]
            );

        if (master_)
        {
            if (shadowPatchField.master() == true)
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
            if (shadowPatchField.master() == false)
            {
                FatalErrorIn
                (
                    "void solidContactFvPatchScalarField::"
                    "calcZoneToZones() const"
                )   << "There is no master patch!" << abort(FatalError);
            }
        }

        if (master_)
        {
            // Create interpolation for patches
            zoneToZones_.set
            (
                shadPatchI,
                new newGgiStandAlonePatchInterpolation
                (
                    zone().globalPatch(),
                    shadowZones()[shadPatchI].globalPatch(),
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
                "void solid4ContactFvPatchVectorField::"
                "calcZoneToZones() const"
            )   << "Attempting to create GGIInterpolation on a shadow patch"
                << abort(FatalError);
        }
    }
}


void Foam::solid4ContactFvPatchVectorField::calcContactPerShadow() const
{
    if (contactPerShadow_.size() > 0)
    {
        FatalErrorIn
        (
            "void thermalContactFvPatchScalarField::"
            "calcContactPerShadow() const"
        )   << "already calculated"
            << abort(FatalError);
    }

    contactPerShadow_.setSize(shadowPatchNames().size());

    forAll(contactPerShadow_, i)
    {
        contactPerShadow_.set
        (
            i,
            new scalarField(patch().size(), 0.0)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::globalPolyPatch&
Foam::solid4ContactFvPatchVectorField::zone() const
{
    if (master_)
    {
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        const solid4ContactFvPatchVectorField& shadowPatchField =
            refCast<const solid4ContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );

        return shadowPatchField.zone();
    }
}


Foam::globalPolyPatch& Foam::solid4ContactFvPatchVectorField::zone()
{
    if (master_)
    {
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        solid4ContactFvPatchVectorField& shadowPatchField =
            const_cast<solid4ContactFvPatchVectorField&>
            (
                refCast<const solid4ContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.zone();
    }
}


const Foam::PtrList<Foam::globalPolyPatch>&
Foam::solid4ContactFvPatchVectorField::shadowZones() const
{
    if (master_)
    {
        if (shadowZones_.empty())
        {
            calcShadowZones();
        }

        return shadowZones_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        const solid4ContactFvPatchVectorField& shadowPatchField =
            refCast<const solid4ContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );

        return shadowPatchField.shadowZones();
    }
}


Foam::PtrList<Foam::globalPolyPatch>&
Foam::solid4ContactFvPatchVectorField::shadowZones()
{
    if (master_)
    {
        if (shadowZones_.empty())
        {
            calcShadowZones();
        }

        return shadowZones_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        solid4ContactFvPatchVectorField& shadowPatchField =
            const_cast<solid4ContactFvPatchVectorField&>
            (
                refCast<const solid4ContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.shadowZones();
    }
}


const Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
Foam::solid4ContactFvPatchVectorField::zoneToZones() const
{
    if (master_)
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

        const solid4ContactFvPatchVectorField& shadowPatchField =
            refCast<const solid4ContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );

        return shadowPatchField.zoneToZones();
    }
}


Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
Foam::solid4ContactFvPatchVectorField::zoneToZones()
{
    if (master_)
    {
        if (zoneToZones_.empty())
        {
            calcZoneToZones();
        }

        return zoneToZones_;
    }
    else
    {
        // We will const_cast the shadow patch so we can delete the weights when
        // the zones move
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        solid4ContactFvPatchVectorField& shadowPatchField =
            const_cast<solid4ContactFvPatchVectorField&>
            (
                refCast<const solid4ContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.zoneToZones();
    }

}


const Foam::newGgiStandAlonePatchInterpolation&
Foam::solid4ContactFvPatchVectorField::zoneToZoneForThisSlave() const
{
    if (master_)
    {
        FatalErrorIn
        (
            "const Foam::newGgiStandAlonePatchInterpolation&"
            " Foam::solid4ContactFvPatchVectorField::"
            "zoneToZoneForThisSlave() const"
        )   << "The master patch is not allowed to call this function"
            << abort(FatalError);
    }

    // The master may have multiple slaves so we need to find which zoneToZone
    // corresponds to the current slave patch
    const wordList& shadPatchNames = shadowPatchField().shadowPatchNames();
    label masterShadowID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterShadowID = shadPatchI;
            break;
        }
    }

    if (masterShadowID == -1)
    {
        FatalErrorIn
        (
            "const Foam::newGgiStandAlonePatchInterpolation&"
            " Foam::solid4ContactFvPatchVectorField::"
            "zoneToZoneForThisSlave() const"
        )   << "Something went wrong when looking for the shadowPatch"
            << abort(FatalError);
    }

    // Return the zoneToZone between the master and the current patch
    return zoneToZones()[masterShadowID];
}


const Foam::globalPolyPatch&
Foam::solid4ContactFvPatchVectorField::zoneForThisSlave() const
{
    if (master_)
    {
        FatalErrorIn
        (
            "const Foam::globalPolyPatch&"
            " Foam::solid4ContactFvPatchVectorField::zoneForThisSlave() const"
        )   << "The master patch is not allowed to call this function"
            << abort(FatalError);
    }

    // The master may have multiple slaves so we need to find which shadowZone
    // corresponds to the current slave patch
    const wordList& shadPatchNames = shadowPatchField().shadowPatchNames();
    label masterShadowID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterShadowID = shadPatchI;
            break;
        }
    }

    if (masterShadowID == -1)
    {
        FatalErrorIn
        (
            "const Foam::globalPolyPatch&"
            " Foam::solid4ContactFvPatchVectorField::zoneForThisSlave() const"
        )   << "Something went wrong when looking for the shadowPatch"
            << abort(FatalError);
    }

    // Return the zoneToZone between the master and the current patch
    return shadowZones()[masterShadowID];
}


// ************************************************************************* //
