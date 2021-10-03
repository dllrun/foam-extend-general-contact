/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "thermalGeneralContactFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


void Foam::thermalGeneralContactFvPatchScalarField::checkConsistentMaster() const
{
	Info<<"In thermalGeneralContact::checkConsistentMaster() line:"<<__LINE__<<endl;
	Info<<"master_ in thermalGeneralContact::checkConsistentMaster(): "<<master_<<endl;
    Info<<"solid4ContactPatch().master() in thermalGeneralContact::checkConsistentMaster(): "<<solid4ContactPatch().master()<<endl;
	if (master_ != solid4ContactPatch().master())
    {
        FatalErrorIn("void checkConsistentMaster() const")
            << "The solid4Contact master patch should be the same as the "
            << "thermalGeneralContact master patch!" << nl
            << abort(FatalError);
    }
	Info<<"In thermalGeneralContact::checkConsistentMaster() line:"<<__LINE__<<endl;
}


const Foam::thermalGeneralContactFvPatchScalarField&
Foam::thermalGeneralContactFvPatchScalarField::shadowPatchField() const
{
    const labelList& shadowPatchIndices =
        solid4ContactPatch().shadowPatchIndices();

    if (shadowPatchIndices.size() != 1)
    {
        FatalErrorIn
        (
            "const Foam::thermalGeneralContactFvPatchScalarField&\n"
            "Foam::thermalGeneralContactFvPatchScalarField::shadowPatchField() const"
        )   << "This function can only be called for a patch with 1 shadow "
            << "patch; this patch has " << shadowPatchIndices.size()
            << " shadow patches!" << abort(FatalError);
    }

    return shadowPatchField(0);
}


const Foam::thermalGeneralContactFvPatchScalarField&
Foam::thermalGeneralContactFvPatchScalarField::shadowPatchField
(
    const label shadI
) const
{
    if (shadI < 0)
    {
        FatalErrorIn
        (
            "const Foam::thermalGeneralContactFvPatchScalarField&\n"
            "Foam::thermalGeneralContactFvPatchScalarField::"
            "shadowPatchField(const label shadI) const"
        )   << "shadI must be a non-negative number!"
            << abort(FatalError);
    }

    const labelList& shadowPatchIndices =
        solid4ContactPatch().shadowPatchIndices();

    if (shadI >= shadowPatchIndices.size())
    {
        FatalErrorIn
        (
            "const Foam::thermalGeneralContactFvPatchScalarField&\n"
            "Foam::thermalGeneralContactFvPatchScalarField::"
            "shadowPatchField(const label shadI) const"
        )   << "shadI is " << shadI << " but this patch has only "
            << shadowPatchIndices.size() << " shadow patches!"
            << abort(FatalError);
    }

    const volScalarField& field =
        db().lookupObject<volScalarField>(dimensionedInternalField().name());

    return
        refCast<const thermalGeneralContactFvPatchScalarField>
        (
            field.boundaryField()[shadowPatchIndices[shadI]]
        );
}


const Foam::solid4ContactFvPatchVectorField&
Foam::thermalGeneralContactFvPatchScalarField::solid4ContactPatch() const
{
	Info<<"In thermalGeneralContact::solid4ContactPatch() line:"<<__LINE__<<endl;
    return
        refCast<const solid4ContactFvPatchVectorField>
        (
            patch().lookupPatchField<volVectorField, vector>(DUName_)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalGeneralContactFvPatchScalarField::thermalGeneralContactFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    master_(false),
    dict_(),
    underRelaxation_(1),
    alpha_(p.size(), 0),
    Tinf_(0),
    Rc_(0),
    beta_(0),
    UTS_(0),
    DUName_("undefined"),
    curTimeIndex_(-1)
{
	Info<<"In thermalGeneralContact c1(p,iF) line:"<<__LINE__<<endl;
    fvPatchScalarField::operator=(patchInternalField());
    gradient() = 0.0;
}


Foam::thermalGeneralContactFvPatchScalarField::thermalGeneralContactFvPatchScalarField
(
    const thermalGeneralContactFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    master_(ptf.master_),
    dict_(ptf.dict_),
    underRelaxation_(ptf.underRelaxation_),
    alpha_(ptf.alpha_, mapper),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_)
{
	Info<<"In thermalGeneralContact c2(ptf,p,iF,mapper) line:"<<__LINE__<<endl;
}


Foam::thermalGeneralContactFvPatchScalarField::thermalGeneralContactFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    master_(dict.lookupOrDefault<Switch>("master", false)),
    dict_(dict),
    underRelaxation_(1),
    alpha_("alpha", dict, p.size()),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    Rc_(0.0),
    beta_(0.0),
    UTS_(0.0),
    DUName_(dict.lookupOrDefault<word>("DUName", "U")), //("DUName", "DU"))
    curTimeIndex_(-1)
{
	Info<<"In thermalGeneralContact c3(p,iF,dict) line:"<<__LINE__<<endl;
    Info<<"master() in thermalGeneralContact c3(p,iF,dict): "<<master()<<endl;
	if (debug)
    {
        Info<< patch().name() << ": " << type() << endl;
    }

    // Read only on master
    if (master())
    {
        underRelaxation_ = readScalar(dict.lookup("underRelaxation"));

        Rc_ = readScalar(dict.lookup("Rc"));

        beta_ = dict.lookupOrDefault<scalar>("beta", 0.1);

        UTS_ = dict.lookupOrDefault<scalar>("UTS", 2e9);

        if (Rc_ < SMALL)
        {
            WarningIn(type() + "::" +  type() + "(...)")
                << "Contact conductivity resistance cannot be exactly zero"
                << nl << "Rc = " << SMALL << " is assumed" << endl;

            Rc_ = SMALL;
        }
    }

    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        gradient() = 0.0;
    }

    if (dict.found("value"))
    {
        Field<scalar>::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        Field<scalar>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }
}


Foam::thermalGeneralContactFvPatchScalarField::thermalGeneralContactFvPatchScalarField
(
    const thermalGeneralContactFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    master_(ptf.master_),
    dict_(ptf.dict_),
    underRelaxation_(ptf.underRelaxation_),
    alpha_(ptf.alpha_),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_)
{
	Info<<"In thermalGeneralContact c4(ptf) line:"<<__LINE__<<endl;
}


Foam::thermalGeneralContactFvPatchScalarField::thermalGeneralContactFvPatchScalarField
(
    const thermalGeneralContactFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    master_(ptf.master_),
    dict_(ptf.dict_),
    underRelaxation_(ptf.underRelaxation_),
    alpha_(ptf.alpha_),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_)
{
	Info<<"In thermalGeneralContact c5(ptf,iF) line:"<<__LINE__<<endl;
}



// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::thermalGeneralContactFvPatchScalarField::~thermalGeneralContactFvPatchScalarField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::wordList&
Foam::thermalGeneralContactFvPatchScalarField::shadowPatchNames() const
{
    return solid4ContactPatch().shadowPatchNames();
}


const Foam::labelList&
Foam::thermalGeneralContactFvPatchScalarField::shadowPatchIndices() const
{
    return solid4ContactPatch().shadowPatchIndices();
}


const Foam::scalarField& Foam::thermalGeneralContactFvPatchScalarField::contact() const
{
    return solid4ContactPatch().contact();
}


Foam::scalar Foam::thermalGeneralContactFvPatchScalarField::underRelaxation() const
{
    if (master())
    {
        return underRelaxation_;
    }
    else
    {
        return shadowPatchField().underRelaxation();
    }
}


Foam::scalar Foam::thermalGeneralContactFvPatchScalarField::Rc() const
{
    if (master())
    {
        return Rc_;
    }
    else
    {
        return shadowPatchField().Rc();
    }
}


Foam::scalar Foam::thermalGeneralContactFvPatchScalarField::beta() const
{
    if (master())
    {
        return beta_;
    }
    else
    {
        return shadowPatchField().beta();
    }
}


Foam::scalar Foam::thermalGeneralContactFvPatchScalarField::UTS() const
{
    if (master())
    {
        return UTS_;
    }
    else
    {
        return shadowPatchField().UTS();
    }
}


Foam::tmp<Foam::scalarField> Foam::thermalGeneralContactFvPatchScalarField::Hc() const
{
    checkConsistentMaster();

    tmp<scalarField> tHc
    (
        new scalarField(patch().size(), 0.0)
    );

    scalarField& Hc = tHc();

    // Contact resistance dependent on contact pressure

    // Pressure dependent contact conductance, where
    // h = hRef*((p/H)**beta)
    // where
    // p = contact pressure
    // H = a measure of the hardness of softer material in Pascal
    // hRef, beta = experimentally fit coefficients
    // beta determines pressure sensitivity
    // beta > 1 very sensitive
    // beta < 0.01 very insensitive
    // See reference: P. Wriggers, Computational Contact Mechanics, Wiley, 2002.

    // Contact resistance is the reciprocal of contact conductance

    // Lookup contact field from DU solid4Contact patch

    // Unit normals
    // Note: the mesh should be in the deformed position so these should be the
    // deformed configuration normals
    const vectorField n = patch().nf();

    // Calculate contact pressure and limit to avoid division by zero in pow
    const scalarField contactPressure =
        max(-n & solid4ContactPatch().traction(), SMALL);

    // Hmnn formula says use Vicker's hardness, but surely it should be
    // Vicker's hardness by 1e6
    // as Vicker's hardness = 0.3*UTS in MPa

    Hc = Foam::pow(contactPressure/(0.3*UTS()), beta())/Rc();

    return tHc;
}


Foam::tmp<Foam::scalarField>
Foam::thermalGeneralContactFvPatchScalarField::frictionFluxRateForThisPatch() const
{
    checkConsistentMaster();

    tmp<scalarField> tfrictionFluxRateForThisPatch
    (
        new scalarField(patch().size(), 0.0)
    );

    // Check the DU patch is of type solid4Contact
    if
    (
        solid4ContactPatch().type()
     != solid4ContactFvPatchVectorField::typeName
    )
    {
        FatalErrorIn
        (
            "thermalGeneralContactFvPatchScalarField::frictionFluxRateForThisPatch()"
        )   << "DU patch " << patch().boundaryMesh()[patch().index()].name()
            << " should be of type solid4Contact "
            << abort(FatalError);
    }
    else
    {
        // Heat flux generated due to friction, stored on the master surface
        scalarField fricFlux(patch().size(), 0.0);

        if (master_)
        {
            fricFlux = solid4ContactPatch().frictionHeatRate();
        }
        else
        {
            // Heat flux on the master global patch
            const scalarField masterZoneFricFlux =
                solid4ContactPatch().shadowPatchField().zone().patchFaceToGlobal
                (
                    solid4ContactPatch().shadowPatchField().frictionHeatRate()()
                );

            // Interpolate the heat flux from the master global patch to the
            // shadow global patch
            const scalarField shadowZoneFricFlux =
                solid4ContactPatch().zoneToZoneForThisSlave().masterToSlave
                (
                    masterZoneFricFlux
                );

            // Convert the global shadow patch to the shadow patch
            // Note: solid4ContactPatch().zone() always returns the master zone
            // and solid4ContactPatch().shadowZones() always returns the slave
            // zones
            fricFlux =
                solid4ContactPatch().zoneForThisSlave().globalFaceToPatch
                (
                    shadowZoneFricFlux
                );
        }

        // In the contact region:
        // masterFlux + slaveFlux - fricFlux = 0
        // And
        // masterFlux = hBar*(slaveTemp - masterTemp) + masterFricFlux
        // Where
        // masterFricFlux = fricFlux*(masterK/masterKsi)
        // masterKsi =
        //     sqrt(slaveRhoC/masterRhoC)*sqrt(masterK*slaveK) + masterK
        // K = conductivity
        // rho = density
        // C = specific heat capacity
        // fricFlux = lossCoeff*shearTraction*relativeVelocity
        // Here we assume the lossCoeff is zero

        // Lookup the thermal conductivity field
		//Note: k is replaced by DT (e.g. you have DT in createFields.H of elasticThermalSolid)
        const volScalarField& k =
            db().lookupObject<volScalarField>("DT");

        // Lookup density times specific heat field
        const volScalarField& rhoC =
            db().lookupObject<volScalarField>("(rho*C)");

        // Current patch index
        const label currentPatchID = patch().index();

        // Get k on the current patch
        const scalarField& curPatchK = k.boundaryField()[currentPatchID];

        // Get rhoC on the current patch
        const scalarField& curPatchRhoC = rhoC.boundaryField()[currentPatchID];

        // For the slave, shadowPatchNames will be length 1, whereas it could be
        // longer for the master
        const wordList shadowPatchNames = this->shadowPatchNames();

        // Reset to zero as we will accumulate frictionFluxRateForThisPatch for
        // all shadow patch
        tfrictionFluxRateForThisPatch() = 0.0;

        forAll(shadowPatchNames, shadI)
        {
            // Get k on the shadow patch
            const scalarField& shadowPatchK =
                k.boundaryField()[shadowPatchIndices()[shadI]];

            // Get rhoC on the shadow patch
            const scalarField& shadowPatchRhoC =
                rhoC.boundaryField()[shadowPatchIndices()[shadI]];

            // Interpolate shadow fields to the curPatch

            // Note: solid4ContactPatch().zone() always returns the master zone
            // and solid4ContactPatch().shadowZones() always returns the slave
            // zones

            scalarField shadowPatchKOnCurPatch;
            scalarField shadowPatchRhoCOnCurPatch;

            if (master_)
            {
                const scalarField shadowZoneK =
                    solid4ContactPatch().shadowZones()[shadI].patchFaceToGlobal
                    (
                        shadowPatchK
                    );

                const scalarField shadowZoneRhoC =
                    solid4ContactPatch().shadowZones()[shadI].patchFaceToGlobal
                    (
                        shadowPatchRhoC
                    );

                const scalarField shadowZoneKOnCurPatch =
                    solid4ContactPatch().zoneToZones()[shadI].slaveToMaster
                    (
                        shadowZoneK
                    );

                const scalarField shadowZoneRhoCOnCurPatch =
                    solid4ContactPatch().zoneToZones()[shadI].slaveToMaster
                    (
                        shadowZoneRhoC
                    );
                shadowPatchKOnCurPatch =
                    solid4ContactPatch().zone().globalFaceToPatch
                    (
                        shadowZoneKOnCurPatch
                    );

                shadowPatchRhoCOnCurPatch =
                    solid4ContactPatch().zone().globalFaceToPatch
                    (
                        shadowZoneRhoCOnCurPatch
                    );
            }
            else
            {
                const scalarField shadowZoneK =
                    solid4ContactPatch().zone().patchFaceToGlobal
                    (
                        shadowPatchK
                    );

                const scalarField shadowZoneRhoC =
                    solid4ContactPatch().zone().patchFaceToGlobal
                    (
                        shadowPatchRhoC
                    );

                const scalarField shadowZoneKOnCurPatch =
                    solid4ContactPatch().zoneToZoneForThisSlave().masterToSlave
                    (
                        shadowZoneK
                    );

                const scalarField shadowZoneRhoCOnCurPatch =
                    solid4ContactPatch().zoneToZoneForThisSlave().masterToSlave
                    (
                        shadowZoneRhoC
                    );

                shadowPatchKOnCurPatch =
                    solid4ContactPatch().shadowZones()[shadI].globalFaceToPatch
                    (
                        shadowZoneKOnCurPatch
                    );

                shadowPatchRhoCOnCurPatch =
                    solid4ContactPatch().shadowZones()[shadI].globalFaceToPatch
                    (
                        shadowZoneRhoCOnCurPatch
                    );
            }

            const scalarField curPatchKsi =
                Foam::sqrt((shadowPatchRhoCOnCurPatch/curPatchRhoC)
                *curPatchK*shadowPatchKOnCurPatch)
              + curPatchK;

            tfrictionFluxRateForThisPatch() +=
                solid4ContactPatch().contactPerShadow()[shadI]
                *fricFlux*(curPatchK/curPatchKsi);
        }
    }

    return tfrictionFluxRateForThisPatch;
}


void Foam::thermalGeneralContactFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);

    alpha_.autoMap(m);
}


void Foam::thermalGeneralContactFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const thermalGeneralContactFvPatchScalarField& dmptf =
        refCast<const thermalGeneralContactFvPatchScalarField>(ptf);

    alpha_.rmap(dmptf.alpha_, addr);
}


void Foam::thermalGeneralContactFvPatchScalarField::updateCoeffs()
{
	Info<<"In thermalGeneralContact::updateCoeffs() line:"<<__LINE__<<endl;
    if (updated())
    {
        return;
    }
	
    // Delete ggi interpolation at the begining of each time step
    if (curTimeIndex_ != db().time().timeIndex())
    {
		Info<<"In thermalGeneralContact::updateCoeffs() line:"<<__LINE__<<endl;
        checkConsistentMaster();
		Info<<"In thermalGeneralContact::updateCoeffs() line:"<<__LINE__<<endl;
        curTimeIndex_ = db().time().timeIndex();
    }

    // Lookup the T field
    // const volScalarField& field =
    //     db().lookupObject<volScalarField>
    //     (
    //         this->dimensionedInternalField().name()
    //     );
	
	Info<<"In thermalGeneralContact::updateCoeffs() line:"<<__LINE__<<endl;
    // Lookup the thermal conductivity field
	//Note: k is replaced by DT (e.g. you have DT in createFields.H of elasticThermalSolid)
    const volScalarField& k = db().lookupObject<volScalarField>("DT"); 
	Info<<"In thermalGeneralContact::updateCoeffs() line:"<<__LINE__<<endl;
	
    // K on the current patch
    const scalarField& curPatchK = k.boundaryField()[patch().index()];

    // Loop through all shadow patches
    // Note: a slave patch will have only one shadow patch (i.e. the master),
    // whereas the master can have multiple shadow patches (i.e. multiple
    // slaves)
    const wordList& shadowPatchNames = this->shadowPatchNames();

    // Accumulate the normal gradient field in the contact areas
    scalarField curPatchSnGradInContactArea(patch().size(), 0.0);

    forAll(shadowPatchNames, shadI)
    {
        // Create the shadow zone temperature field
        scalarField shadowPatchTOnCurPatch;
        if (master())
        {
            // Note: solid4ContactPatch().zone() always returns the master zone
            // and solid4ContactPatch().shadowZones() always returns the slave
            // zones
            const scalarField shadowZoneT =
                solid4ContactPatch().shadowZones()[shadI].patchFaceToGlobal
                (
                    shadowPatchField(shadI)
                );

            // Interpolate shadow temperature field to the current zone
            const scalarField shadowZoneTOnCurPatch =
                solid4ContactPatch().zoneToZones()[shadI].slaveToMaster
                (
                    shadowZoneT
                );

            // Create the temperature field for this patch
            shadowPatchTOnCurPatch =
                solid4ContactPatch().zone().globalFaceToPatch
                (
                    shadowZoneTOnCurPatch
                );
        }
        else
        {
            // Note: zone() is the master global patch
            // Note 2: the slave always has only one master i.e.
            // one shadowPatchField
            const scalarField shadowZoneT =
                solid4ContactPatch().zone().patchFaceToGlobal
                (
                    shadowPatchField()
                );

            // Interpolate shadow temperature field to the current zone
            const scalarField shadowZoneTOnCurPatch =
                solid4ContactPatch().zoneToZoneForThisSlave().masterToSlave
                (
                    shadowZoneT
                );

            // Create the temperature field for this patch
            shadowPatchTOnCurPatch =
                solid4ContactPatch().shadowZones()[shadI].globalFaceToPatch
                (
                    shadowZoneTOnCurPatch
                );
        }

        // Calculate current contact conductance
        const scalarField curPatchH = Hc();

        // Calculate the heat flux through the current patch (in the contact
        // area)
        const scalarField curPatchFluxInContactArea =
          - curPatchH*(shadowPatchTOnCurPatch - *this)
          - frictionFluxRateForThisPatch();

        // Get the contact indicator field for this contact pair
        const scalarField contactPerShadow =
            solid4ContactPatch().contactPerShadow()[shadI];

        // Convert flux to normal gradient
        // fluxInContactArea == k*snGradTInContactArea
        // therefore
        // snGradTInContactArea == fluxInContactArea/k

        // Accumulate for all contact pairs
        curPatchSnGradInContactArea +=
            -contactPerShadow*curPatchFluxInContactArea/curPatchK;
    }

    // Next, we will add the thermal convections contributions for regions not
    // in contact:
    // k*snGradTNotInContact = -alpha*(T - Tinf)
    // i.e. heat flux within solid == heat flux due to convection at the
    // surface, therefore:
    // snGradTsnGradTNotInContact = -(alpha/k)*(T - Tinf)
    // Info<< nl
    //     << solid4ContactPatch().contact().size() << nl
    //     << alpha_.size() << nl
    //     << this->size() << nl
    //     << curPatchK.size() << endl;
    const scalarField curPatchPatchSnGradNotInContact =
        (1.0 - solid4ContactPatch().contact())
       *(-alpha_/curPatchK)*(*this - Tinf_);

    // Merge contributions for "in contact" and "not in contact" regions
    const scalarField curPatchPatchSnGrad =
        curPatchSnGradInContactArea + curPatchPatchSnGradNotInContact;

    // Set gradient using under-relaxation
    gradient() =
        underRelaxation()*curPatchPatchSnGrad
      + (1.0 - underRelaxation())*gradient();

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::thermalGeneralContactFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
	Info<<"In thermalGeneralContact::evaluate(..) line:"<<__LINE__<<endl;
    if (!this->updated())
    {
        this->updateCoeffs();
    }

	Info<<"In thermalGeneralContact::evaluate(..) line:"<<__LINE__<<endl;
    const fvPatchField<vector>& gradT =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );
	
	Info<<"In thermalGeneralContact::evaluate(..) line:"<<__LINE__<<endl;
	
    // Unit normals
    const vectorField n = patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // Correction vectors
    const vectorField k = (I - sqr(n)) & delta;

    Field<scalar>::operator=
    (
        patchInternalField()
      + (k & gradT.patchInternalField())
      + gradient()/patch().deltaCoeffs()
    );

    fvPatchField<scalar>::evaluate();
}


void Foam::thermalGeneralContactFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;
    alpha_.writeEntry("alpha", os);
    os.writeKeyword("Tinf")
        << Tinf_ << token::END_STATEMENT << nl;

    if (master())
    {
        os.writeKeyword("underRelaxation") << underRelaxation_
            << token::END_STATEMENT << nl;
        os.writeKeyword("Rc")
            << Rc_ << token::END_STATEMENT << nl;
        os.writeKeyword("beta")
            << beta_ << token::END_STATEMENT << nl;
        os.writeKeyword("UTS")
            << UTS_ << token::END_STATEMENT << nl;
    }

    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        thermalGeneralContactFvPatchScalarField
    );
}

// ************************************************************************* //
