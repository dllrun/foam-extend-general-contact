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

InClass
    solid4GeneralContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solid4GeneralContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
//#include "lookupSolidModel.H"
//#define __LINE__ 0
#define ISDEBUG false
#define ActivePairDEBUG false
#define currentMasterDEBUG false
#define localSlaveDEBUG false
#define normalModelDEBUG false
#define masterOfPairDEBUG false
#define bbOffDEBUG false


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// *********************** START solid GeneralContact ************************
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
// *********************** END solid GeneralContact ************************

/*
bool Foam::solid4GeneralContactFvPatchVectorField::movingMesh() const
{
	
    // Check if the solid model moves the mesh
    return lookupSolidModel(patch().boundaryMesh().mesh()).movingMesh();
	
}
*/

//******************* based on solid General******************
bool Foam::solid4GeneralContactFvPatchVectorField::firstPatchInList() const
{
    if (!firstPatchPtr_)
    {
        calcFirstPatchInList();
    }
	/*
	Info<< "The current field in firstPatchInList() is "<< dimensionedInternalField().name()<< endl;
	Info<< "The current patch in firstPatchInList()) is "<< patch().name()<< endl;			
	Info<< "The size of current patch in firstPatchInList() is "<< patch().size()<< endl;
	Info<< "*firstPatchPtr_ in firstPatchInList() is "<<*firstPatchPtr_<< endl;
    */
	return *firstPatchPtr_;
}


void Foam::solid4GeneralContactFvPatchVectorField::calcFirstPatchInList() const   //// CHECK method to calculate global master
{
	
    if (firstPatchPtr_)
    {
        FatalErrorIn
            (
                "void Foam::solid4GeneralContactFvPatchVectorField::"
                "calcFirstPatchInList() const"
            )   << "firstPatchPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solid4GeneralContact patch i.e. the one
    // with the lowest patch index
	/*
	Info<< "In calcFirstPatchInList() "<<__LINE__<<endl;
	Info<< "firstPatchIndexInList() in calcFirstPatchInList() "<<firstPatchIndexInList()<<endl;
	Info<< "patch().index() in calcFirstPatchInList() "<<patch().index()<<endl;
    Info<< "patch().name() in calcFirstPatchInList() "<<patch().name()<<endl;
	*/
	if (firstPatchIndexInList() == patch().index())
    {
        firstPatchPtr_ = new bool(true);
    }
    else
    {
        firstPatchPtr_ = new bool(false);
    }
	
}

Foam::label Foam::solid4GeneralContactFvPatchVectorField::firstPatchIndexInList()
const
{
    if (!firstPatchIndexInListPtr_)
    {
        calcFirstPatchIndexInList();
    }


    return *firstPatchIndexInListPtr_;
}

void Foam::solid4GeneralContactFvPatchVectorField::calcFirstPatchIndexInList() const
{
    if (firstPatchIndexInListPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solid4GeneralContactFvPatchVectorField::"
            "calcFirstPatchIndexInList() const"
        )   << "firstPatchIndexInListPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solid4GeneralContact patch i.e. the one
    // with the lowest patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    firstPatchIndexInListPtr_ = new label(-1);
    label& gMasterID = *firstPatchIndexInListPtr_;

    forAll(field.boundaryField(), patchI)
    {
		Info<<"In calcGlobalMaster line: "<<__LINE__<<endl;
        if
        (
            field.boundaryField()[patchI].type()
            == solid4GeneralContactFvPatchVectorField::typeName
        )
        {
			Info<<"In calcGlobalMaster line: "<<__LINE__<<endl;
            gMasterID = patchI;

            break;
        }		
    }
	Info<<"In calcGlobalMaster line: "<<__LINE__<<endl;
		Info<<"gMasterID in calcGlobalMaster: "<<gMasterID<<endl;

    // Check there is only one global master

    label GMasterID = returnReduce(gMasterID, maxOp<label>());

    if (gMasterID != GMasterID)
    {
        FatalErrorIn
        (
            "solid4GeneralContactFvPatchVectorField::"
            "calcFirstPatchIndexInList() const"
        )   << "There are multiple global masters" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< nl << "The global master contact patch is "
            << patch().boundaryMesh()[gMasterID].name() << endl;
    }
}

//******************* current master **********************
bool Foam::solid4GeneralContactFvPatchVectorField::currentMaster() const
{
    if (!currentMasterPtr_)
    {
        calcCurrentMaster();
    }
	#if(currentMasterDEBUG) 
	Info<< "In currentMaster() "<<__LINE__<<endl;
	Info<< "patch().index() in currentMaster() "<<patch().index()<<endl;
	Info<< "patch().name() in currentMaster() "<<patch().name()<<endl;
	Info<< "*currentMasterPtr_ in currentMaster() "<<*currentMasterPtr_<<endl;
	#endif
	
	return *currentMasterPtr_;
}


void Foam::solid4GeneralContactFvPatchVectorField::calcCurrentMaster() const   //// CHECK method to calculate global master
{
	
    if (currentMasterPtr_)
    {
        FatalErrorIn
            (
                "void Foam::solid4GeneralContactFvPatchVectorField::"
                "calcCurrentMaster() const"
            )   << "currentMasterPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solid4GeneralContact patch i.e. the one
    // with the lowest patch index
	
	//const boolList& locSlave = localSlave();
	const boolList& locSlave = localTSlave();
	#if(currentMasterDEBUG) 
	Info<< "In calcCurrentMaster() "<<__LINE__<<endl;
	Info<< "patch().index() in calcCurrentMaster() "<<patch().index()<<endl;
    #endif
	
	forAll(locSlave, slaveI)
    {
        if (locSlave[slaveI])
        {
			currentMasterPtr_ = new bool(true);
		}
		else
		{
			currentMasterPtr_ = new bool(false);
		}
    }	
		
}

//***************** end current master ********************


//******************* current slave **********************
bool Foam::solid4GeneralContactFvPatchVectorField::currentSlave() const
{
    if (!currentSlavePtr_)
    {
        calcCurrentSlave();
    }
	Info<< "In currentSlave() "<<__LINE__<<endl;
	Info<< "patch().index() in currentSlave() "<<patch().index()<<endl;
	Info<< "patch().name() in currentSlave() "<<patch().name()<<endl;
	return *currentSlavePtr_;
}


void Foam::solid4GeneralContactFvPatchVectorField::calcCurrentSlave() const   //// CHECK method to calculate global master
{
	
    if (currentSlavePtr_)
    {
        FatalErrorIn
            (
                "void Foam::solid4GeneralContactFvPatchVectorField::"
                "calcCurrentSlave() const"
            )   << "currentSlavePtr_ already set" << abort(FatalError);
    }

    // The global master is the first solid4GeneralContact patch i.e. the one
    // with the lowest patch index 
	
	//const List<Foam::word>& curSlave = slavePatchNames();
	
	Info<< "In calcCurrentSlave() "<<__LINE__<<endl;
	Info<< "patch().index() in calcCurrentSlave() "<<patch().index()<<endl;
    /*
	forAll(curSlave, slaveI)
    {
        if (curSlave[slaveI])
        {
			currentSlavePtr_ = new bool(true);
		}
		else
		{
			currentSlavePtr_ = new bool(false);
		}
    }
	*/	
		
}

//***************** end current slave ********************



const Foam::boolList&
Foam::solid4GeneralContactFvPatchVectorField::localSlave() const
{
    if (!localSlavePtr_)
    {
        calcLocalSlave();
    }
	
	#if(localSlaveDEBUG)
	Info<< "The current field in localSlave() is "<< dimensionedInternalField().name()<< endl;
	Info<< "The current patch in localSlave() is "<< patch().name()<< endl;			
	Info<< "The size of current patch in localSlave() is "<< patch().size()<< endl;
	Info<< "*localSlavePtr_ in localSlave() is "<<*localSlavePtr_<< endl;
    #endif
	
	return *localSlavePtr_;
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

    localSlavePtr_ = new boolList(slavePatchNames().size(), false);

    boolList& localSlave = *localSlavePtr_;
	
	Info<< "In calcLocalSlave() "<<__LINE__<<endl;
	Info<< "patch().index() in calcLocalSlave() "<<patch().index()<<endl;
    forAll(localSlave, slaveI)
    {
        if (patch().index() < slavePatchIndices()[slaveI])
        {
            localSlave[slaveI] = true;

            Info<< "solid4GeneralContact: "
                << patch().name() << " (master)" << " to "
                << slavePatchNames()[slaveI] << " (slave)" << endl;
        }
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
		#if(bbOffDEBUG) 
		Info<<"patch().magSf() in calcBbOffset() is "<<patch().magSf()<< endl;
		Info<<"sqrt(patch().magSf()) in calcBbOffset() is "<<sqrt(patch().magSf())<< endl;
		Info<<"min(sqrt(patch().magSf())) in calcBbOffset() is "<<min(sqrt(patch().magSf()))<< endl;
        #endif
		minDim = min(sqrt(patch().magSf()));
    }
	
	#if(bbOffDEBUG)
	Info<<"bbOffset_ in calcBbOffset() is "<<bbOffset_<< endl;
    #endif
	
	bbOffset_ = 6.5*returnReduce(minDim, minOp<scalar>());
	
	#if(bbOffDEBUG)
	Info<<"bbOffset_ in calcBbOffset() is "<<bbOffset_<< endl;
	#endif
	
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
	Info<<"The bbOffset bbOffset() is "<<bbOffset_<< endl;

    return bbOffset_;
}


//slave bbOff functions

void Foam::solid4GeneralContactFvPatchVectorField::calcSlaveBbOffset() const
{
	Info<<"In calcSlaveBbOffset() line:"<<__LINE__<<endl;
    if (!slaveBbOffsetList_.empty())
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::solid4GeneralContactFvPatchVectorField::"
            "calcSlaveBbOffset() const"
        )   << "already set" << abort(FatalError);
    }

    // We will set the BB offset to five times the average dimension of the
    // smallest face on the zone

    scalarField minDimList(slavePatchNames().size(),GREAT);
	
	Info<<"In calcSlaveBbOffset() line:"<<__LINE__<<endl;
    if (patch().size() > 0)
    {
		//const List<Foam::word>& curSlave = slavePatchNames();
		forAll(minDimList, shadPatchI)
		{
			/*	
			patch().patch().boundaryMesh()
                [
                    slavePatchIndices()[0]
                ].meshPoints();
			*/
		
		minDimList[shadPatchI] = min(sqrt(patch().boundaryMesh()[slavePatchIndices()[shadPatchI]].magSf()));
		
		Info<<"In calcSlaveBbOffset() line:"<<__LINE__<<endl;
		Info<<"minDimList[shadPatchI] calcSlaveBbOffset() :"<<minDimList[shadPatchI]<<endl;
		//Info<<"slavePatchNames()[shadPatchI] calcSlaveBbOffset() :"<<slavePatchNames()[shadPatchI]<<endl;
		}
		
    }	
		//Info<<"slaveBbOffsetList_.size() in calcSlaveBbOffset() :"<<slavePatchNames().size()<<endl;
		scalarField slaveBbOff (slavePatchNames().size(),0.0);
			
		Info<<"In calcSlaveBbOffset() line:"<<__LINE__<<endl;
		
		forAll(slaveBbOff, shadPatchI)
		{
			Info<<"In calcSlaveBbOffset() line:"<<__LINE__<<endl;
		slaveBbOff[shadPatchI] = 6.0*returnReduce(minDimList[shadPatchI], minOp<scalar>());
		}
				
		Info<<"In calcSlaveBbOffset() line:"<<__LINE__<<endl;
		Info<<"slaveBbOffsetList_ in calcSlaveBbOffset() is "<<slaveBbOffsetList_<< endl;
		slaveBbOffsetList_ = slaveBbOff;
		Info<<"In calcSlaveBbOffset() line:"<<__LINE__<<endl;
		Info<<"slaveBbOffsetList_ in calcSlaveBbOffset() is "<<slaveBbOffsetList_<< endl;
		
	
    if (debug)
    {
        Info<< nl << "The bbOffset is " << slaveBbOffsetList_ << endl;
    }
}


Foam::scalarField 
Foam::solid4GeneralContactFvPatchVectorField::slaveBbOffset() const
{
    //if (slaveBbOffsetList_ == 0)
	if(slaveBbOffsetList_.empty())
    {
       calcSlaveBbOffset();
    }
	//Info<<"The slavePatchNames_ in slaveBbOffset() is "<<slavePatchNames_<< endl;
	Info<<"The bbOffset slaveBbOffset() is "<<slaveBbOffsetList_<< endl;

    return slaveBbOffsetList_;
}

//******************** END based on solid General*****************




void Foam::solid4GeneralContactFvPatchVectorField::makeSlavePatchNames()const
/*(
    const dictionary& dict
) const */
{
	//Info<<"In makeSlavePatchNames() line:"<<__LINE__<<endl;
	//********************** based on solid General*************
	if (slavePatchNames_ || slavePatchIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solid4GeneralContactFvPatchVectorField::"
            "calcSlavePatchNames() const"
        )   << "slavePatchNames_ or slavePatchIndices_ already set"
            << abort(FatalError);
    }
	
	// Add each solid4GeneralContact patch in the order of increasing patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );
		
	// Count slave patches

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
	
	
	slavePatchNames_ = new wordList(nShadPatches);
	//Info<<"In makeSlavePatchNames() line:"<<__LINE__<<endl;	
	wordList& slavePatchNames = *slavePatchNames_;
	//Info<<"slavePatchNames in makeSlavePatchNames() "<<slavePatchNames<<endl;
	
	slavePatchIndicesPtr_ = new labelList(nShadPatches);
    labelList& slavePatchIndices = *slavePatchIndicesPtr_;
	
	// Record slave patch names

    label shadPatchI = 0;

    forAll(field.boundaryField(), patchI)
    {
		//Info<<"In makeSlavePatchNames() line:"<<__LINE__<<endl;
		//Info<<"patch().index() in makeSlavePatchNames(): "<<patch().index()<<endl;	
        if
        (
            field.boundaryField()[patchI].type()
            == solid4GeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
			//Info<<"shadPatchI in makeSlavePatchNames(): "<<shadPatchI<<endl;	
            slavePatchNames[shadPatchI] = patch().boundaryMesh()[patchI].name();           
			//Info<<"slavePatchNames[shadPatchI] in makeSlavePatchNames(): "<<slavePatchNames[shadPatchI]<<endl;	
			slavePatchIndices[shadPatchI++] = patchI;
			//Info<<"shadPatchI in makeSlavePatchNames(): "<<shadPatchI<<endl;
			//Info<<"patchI in makeSlavePatchNames(): "<<patchI<<endl;
		}
    }
	
	//******************** END based on solid General*************			
}

//****************** Start Test for shadPatchI dependent function **************
void Foam::solid4GeneralContactFvPatchVectorField::calcNormalContactModels
(
    const dictionary& dict
) const
{
	#if(normalModelDEBUG)
	Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
    #endif
	
	normalModels_.setSize(slavePatchNames().size());
	
	#if(normalModelDEBUG)
	Info<<"normalModels_.size() in makeNormalModels(..): "<<normalModels_.size()<<endl;
    #endif
	
	const boolList& locSlave = localSlave();
	
	forAll (normalModels_, shadPatchI)
    {
		/*
        // Check if only one slave patch is defined
        const dictionary* contactDictPtr = NULL;
        if (normalModels_.size() == 1)
        {
			Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
            if (dict.found("generalNormalContactModel") )
            {				
                contactDictPtr = &dict;
				Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
            }
            else
            {
                contactDictPtr =
                    &dict.subDict
                    (
                        patch().name() + "_to_"
                      + slavePatchNames()[shadPatchI] + "_dict"
                    );
            }
        }
        else
        {
            contactDictPtr =
                &dict.subDict
                (
                    patch().name() + "_to_"
                  + slavePatchNames()[shadPatchI] + "_dict"
                );
        }
        const dictionary& contactDict = *contactDictPtr;
		Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
		*/
		#if(normalModelDEBUG)
		Info<<"locSlave[shadPatchI] in makeNormalModels(): "<<locSlave[shadPatchI]<<endl;
		#endif
		
		// Only the local slave creates the contact model
        if (locSlave[shadPatchI])
        {
			const dictionary* contactDictPtr = NULL;
			#if(normalModelDEBUG)
			Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
            #endif
			
			if (dict.found("generalNormalContactModel") )
            {				
                contactDictPtr = &dict;
				#if(normalModelDEBUG)
				Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
				#endif
				
			}
			
			const dictionary& contactDict = *contactDictPtr;
			
			#if(normalModelDEBUG)
			Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
			#endif
			
			// Create contact model
			normalModels_.set
			(
				shadPatchI,
				generalNormalContactModel::New
				(
					word(contactDict.lookup("generalNormalContactModel")),
					patch().boundaryMesh()[slavePatchIndices()[shadPatchI]],
					contactDict,
					patch().index(),                  // master
					slavePatchIndices()[shadPatchI], // slave
					zone().globalPatch(),
					slaveZones()[shadPatchI].globalPatch()
				)
			);
		}
		
		#if(normalModelDEBUG)
		Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
		#endif
	}

    // Initialise penalty scales to -1
    Info<< "    Initialising stored previous normalPenaltyFactors" << endl;
    normalPenaltyFactors_.setSize(normalModels_.size(), -1);
}


//****************** End Test for shadPatchI dependent function **************


void Foam::solid4GeneralContactFvPatchVectorField::makeNormalModels
(
    const dictionary& dict
) const
{
	#if(normalModelDEBUG)
	Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
    #endif
	
	normalModels_.setSize(slavePatchNames().size());
	
	#if(normalModelDEBUG)
	Info<<"normalModels_.size() in makeNormalModels(..): "<<normalModels_.size()<<endl;
    #endif
	
	const boolList& locSlave = localSlave();
	
	forAll (normalModels_, shadPatchI)
    {
		/*
        // Check if only one slave patch is defined
        const dictionary* contactDictPtr = NULL;
        if (normalModels_.size() == 1)
        {
			Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
            if (dict.found("generalNormalContactModel") )
            {				
                contactDictPtr = &dict;
				Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
            }
            else
            {
                contactDictPtr =
                    &dict.subDict
                    (
                        patch().name() + "_to_"
                      + slavePatchNames()[shadPatchI] + "_dict"
                    );
            }
        }
        else
        {
            contactDictPtr =
                &dict.subDict
                (
                    patch().name() + "_to_"
                  + slavePatchNames()[shadPatchI] + "_dict"
                );
        }
        const dictionary& contactDict = *contactDictPtr;
		Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
		*/
		#if(normalModelDEBUG)
		Info<<"locSlave[shadPatchI] in makeNormalModels(): "<<locSlave[shadPatchI]<<endl;
		#endif
		
		// Only the local slave creates the contact model
        if (locSlave[shadPatchI])
        {
			const dictionary* contactDictPtr = NULL;
			#if(normalModelDEBUG)
			Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
            #endif
			
			if (dict.found("generalNormalContactModel") )
            {				
                contactDictPtr = &dict;
				#if(normalModelDEBUG)
				Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
				#endif
				
			}
			
			const dictionary& contactDict = *contactDictPtr;
			
			#if(normalModelDEBUG)
			Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
			#endif
			
			// Create contact model
			normalModels_.set
			(
				shadPatchI,
				generalNormalContactModel::New
				(
					word(contactDict.lookup("generalNormalContactModel")),
					patch(),
					contactDict,
					patch().index(),                  // master
					slavePatchIndices()[shadPatchI], // slave
					zone().globalPatch(),
					slaveZones()[shadPatchI].globalPatch()
				).ptr()
			);
		}
		
		#if(normalModelDEBUG)
		Info<<"In makeNormalModels(..) line:"<<__LINE__<<endl;
		#endif
	}

    // Initialise penalty scales to -1
    Info<< "    Initialising stored previous normalPenaltyFactors" << endl;
    normalPenaltyFactors_.setSize(normalModels_.size(), -1);
}


void Foam::solid4GeneralContactFvPatchVectorField::makeFrictionModels
(
    const dictionary& dict
) const
{
    frictionModels_.setSize(slavePatchNames().size());
	
	const boolList& locSlave = localSlave();
	
    forAll (frictionModels_, shadPatchI)
    {
		/*
        // Check if only one slave patch is defined
        const dictionary* contactDictPtr = NULL;
        if (frictionModels_.size() == 1)
        {
            if (dict.found("generalFrictionContactModel") )
            {
                contactDictPtr = &dict;
            }
            else
            {
                contactDictPtr =
                    &dict.subDict
                    (
                        patch().name() + "_to_"
                      + slavePatchNames()[shadPatchI] + "_dict"
                    );
            }
        }
        else
        {
            contactDictPtr =
                &dict.subDict
                (
                    patch().name() + "_to_"
                  + slavePatchNames()[shadPatchI] + "_dict"
                );
        }
        const dictionary& contactDict = *contactDictPtr;
		*/
		Info<<"locSlave[shadPatchI] in makeFrictionModels(): "<<locSlave[shadPatchI]<<endl;
		if (locSlave[shadPatchI])
        {
			const dictionary* contactDictPtr = NULL;
			if (dict.found("generalFrictionContactModel") )
            {
                contactDictPtr = &dict;
            }
			
			const dictionary& contactDict = *contactDictPtr;
			// Create contact model
			frictionModels_.set
			(
				shadPatchI,
				generalFrictionContactModel::New
				(
					word(contactDict.lookup("generalFrictionContactModel")),
					patch(),
					contactDict,
					patch().index(),                 // master
					slavePatchIndices()[shadPatchI] // slave
				).ptr()
			);
		}
    }
	Info<<"In makeFrictionModels(..) line:"<<__LINE__<<endl;
}


void Foam::solid4GeneralContactFvPatchVectorField::clearOut()
{
	Info<<"In clearOut() line:"<<__LINE__<<endl;
    if (debug)
    {
        InfoIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::clearOut()"
        )   << patch().name() << " : clearOut" << endl;
    }
	//************ based on solid General*************
	deleteDemandDrivenData(firstPatchPtr_);
    deleteDemandDrivenData(firstPatchIndexInListPtr_);
	deleteDemandDrivenData(currentMasterPtr_);
    deleteDemandDrivenData(currentMasterIndexPtr_);
    deleteDemandDrivenData(localSlavePtr_);
	deleteDemandDrivenData(curPatchTractionPtr_);
	//************ END based on solid General*************
    
	deleteDemandDrivenData(slavePatchIndicesPtr_);
    deleteDemandDrivenData(zonePtr_);
    slaveZones_.clear();
    zoneToZones_.clear();
    scaleTractionFieldPtr_.clear();
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solid4GeneralContactFvPatchVectorField::solid4GeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:	
    solidTractionFvPatchVectorField(p, iF),
	//******** based on solid General***************
    firstPatchPtr_(NULL),
    firstPatchIndexInListPtr_(NULL),
	currentMasterPtr_(NULL),
    currentMasterIndexPtr_(NULL),
	currentSlavePtr_(NULL),
	localSlavePtr_(NULL),
	//******** END based on solid General*************
    dict_(),
	//master_(false),
    writeZoneVTK_(false),
    writePointDistanceFields_(false),
    slavePatchNames_(),
    slavePatchIndicesPtr_(NULL),
    rigidMaster_(false),
    normalModels_(),
    frictionModels_(),
    normalPenaltyFactors_(),
    zonePtr_(NULL),
    slaveZones_(),
    zoneToZones_(),
    quickReject_(Foam::newGgiInterpolation::AABB),
    regionOfInterestTopCorner_(vector::max),
    regionOfInterestBottomCorner_(vector::min),
    regionOfInterest_(vector::min, vector::max),
    contact_(0),
    contactPerSlave_(),
    scaleFaceTractionsNearDownstreamPatch_(false),
    scaleTractionFieldPtr_(),
    curTimeIndex_(-1),
	//******** based on solid General*************
	curPatchTractionPtr_(NULL),
	bbOffset_(0.0),
	slaveBbOffsetList_(0.0)
	//******** END based on solid General*************
{
	Info<<"Does it enter here? in C1(p, iF)"<<__LINE__<<endl;
}


Foam::solid4GeneralContactFvPatchVectorField::solid4GeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
	//******** based on solid General***************
    firstPatchPtr_(NULL),
    firstPatchIndexInListPtr_(NULL),
	currentMasterPtr_(NULL),
    currentMasterIndexPtr_(NULL),
	currentSlavePtr_(NULL),
	localSlavePtr_(NULL),
	//******** END based on solid General*************
    dict_(dict),
    //master_(dict.lookup("master")),
    writeZoneVTK_(dict.lookupOrDefault<Switch>("writeZoneVTK", false)),
    writePointDistanceFields_
    (
        dict.lookupOrDefault<Switch>("writePointDistanceFields", false)
    ),
    slavePatchNames_(),
    slavePatchIndicesPtr_(NULL),
    rigidMaster_(false),
    normalModels_(),
    frictionModels_(),
    normalPenaltyFactors_(),
    zonePtr_(NULL),
    slaveZones_(),
    zoneToZones_(),
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
    contact_(patch().size(), 0.0),
    contactPerSlave_(),
    scaleFaceTractionsNearDownstreamPatch_
    (
        dict.lookupOrDefault<Switch>
        (
            "scaleFaceTractionsNearDownstreamPatch",
            Switch(false)
        )
    ),
    scaleTractionFieldPtr_(),
    curTimeIndex_(-1),
	//******** based on solid General*************
	curPatchTractionPtr_(NULL),
	bbOffset_(0.0),
	slaveBbOffsetList_(0.0)
	//******** END based on solid General*************
{
	Info<<"Does it enter here? in C2(p, iF, dict)"<<__LINE__<<endl;
    if (debug)
    {
        Info<< "Creating " << solid4GeneralContactFvPatchVectorField::typeName
            << " patch" << endl;
    }

    // Master creates contact laws
    //if (master_) 
	if (firstPatchPtr_)//if (ocalSlavePtr_)
    {
        rigidMaster_ = Switch(dict.lookup("rigidMaster"));

        if (debug)
        {
            Info<< "    writePointDistanceFields: " << writePointDistanceFields_
                << endl;
        }

        if (scaleFaceTractionsNearDownstreamPatch_)
        {
            WarningIn(type() + "::" + type())
                << "scaleFaceTractionsNearDownstreamPatch can only be applied on"
                << "the slave patch: this option will be ignored for the master "
                << "patch!" << endl;
        }
    }
	
	if (currentMasterPtr_)//if (ocalSlavePtr_)
    {
        rigidMaster_ = Switch(dict.lookup("rigidMaster"));

        if (debug)
        {
            Info<< "    writePointDistanceFields: " << writePointDistanceFields_
                << endl;
        }

        if (scaleFaceTractionsNearDownstreamPatch_)
        {
            WarningIn(type() + "::" + type())
                << "scaleFaceTractionsNearDownstreamPatch can only be applied on"
                << "the slave patch: this option will be ignored for the master "
                << "patch!" << endl;
        }
    }

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


Foam::solid4GeneralContactFvPatchVectorField::solid4GeneralContactFvPatchVectorField
(
    const solid4GeneralContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper),
    //******** based on solid General***************
    firstPatchPtr_(NULL),
    firstPatchIndexInListPtr_(NULL),
	currentMasterPtr_(NULL),
    currentMasterIndexPtr_(NULL),
	currentSlavePtr_(NULL),
	localSlavePtr_(NULL),
	//******** END based on solid General*************
	dict_(ptf.dict_),
    //master_(ptf.master_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    writePointDistanceFields_(ptf.writePointDistanceFields_),
    slavePatchNames_(ptf.slavePatchNames_),
    slavePatchIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    slaveZones_(),
    zoneToZones_(),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    contact_(ptf.contact_),
    contactPerSlave_(),
    scaleFaceTractionsNearDownstreamPatch_
    (
        ptf.scaleFaceTractionsNearDownstreamPatch_
    ),
    scaleTractionFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_),
	//******** based on solid General*************
	curPatchTractionPtr_(NULL),
	bbOffset_(ptf.bbOffset_),
	slaveBbOffsetList_(ptf.bbOffset_)
	//******** END based on solid General*************
{
    Info<<"Does it enter here? in C3(ptf, p, iF, mapper)"<<__LINE__<<endl;
	// Do not copy pointer objects: they will be re-created.
	
	//******************* based on solid General******************
	// Copy pointer objects
    if (ptf.firstPatchPtr_)
    {
        firstPatchPtr_ = new bool(*ptf.firstPatchPtr_);
    }
	
	if (ptf.firstPatchIndexInListPtr_)
    {
        firstPatchIndexInListPtr_ = new label(*ptf.firstPatchIndexInListPtr_);
    }
	
	if (ptf.currentMasterPtr_)
    {
        currentMasterPtr_ = new bool(*ptf.currentMasterPtr_);
    }
	
	if (ptf.currentMasterIndexPtr_)
    {
        currentMasterIndexPtr_ = new label(*ptf.currentMasterIndexPtr_);
    }
	
	if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }
	
	if (ptf.curPatchTractionPtr_)
    {
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }
	//******************* END based on solid General******************
}


Foam::solid4GeneralContactFvPatchVectorField::solid4GeneralContactFvPatchVectorField
(
    const solid4GeneralContactFvPatchVectorField& ptf
)
:
    solidTractionFvPatchVectorField(ptf),
	//******** based on solid General***************
    firstPatchPtr_(NULL),
    firstPatchIndexInListPtr_(NULL),
	currentMasterPtr_(NULL),
    currentMasterIndexPtr_(NULL),
	currentSlavePtr_(NULL),
	localSlavePtr_(NULL),
	//******** END based on solid General*************
    dict_(ptf.dict_),
    //master_(ptf.master_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    writePointDistanceFields_(ptf.writePointDistanceFields_),
    slavePatchNames_(ptf.slavePatchNames_),
    slavePatchIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    slaveZones_(),
    zoneToZones_(),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    contact_(ptf.contact_),
    contactPerSlave_(),
    scaleFaceTractionsNearDownstreamPatch_
    (
        ptf.scaleFaceTractionsNearDownstreamPatch_
    ),
    scaleTractionFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_),
	//******** based on solid General*************
	curPatchTractionPtr_(NULL),
	bbOffset_(ptf.bbOffset_),
	slaveBbOffsetList_(ptf.bbOffset_)
	//******** END based on solid General*************
{
    Info<<"Does it enter here? in C4(ptf)"<<__LINE__<<endl;
	// Do not copy pointer objects
	
	//******************* based on solid General******************
	// Copy pointer objects
	if (ptf.firstPatchPtr_)
    {
        firstPatchPtr_ = new bool(*ptf.firstPatchPtr_);
    }

    if (ptf.firstPatchIndexInListPtr_)
    {
        firstPatchIndexInListPtr_ = new label(*ptf.firstPatchIndexInListPtr_);
    }
	
	if (ptf.currentMasterPtr_)
    {
        currentMasterPtr_ = new bool(*ptf.currentMasterPtr_);
    }

    if (ptf.currentMasterIndexPtr_)
    {
        currentMasterIndexPtr_ = new label(*ptf.currentMasterIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }
	
	if (ptf.curPatchTractionPtr_)
    {
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }
	//******************* END based on solid General******************
}


Foam::solid4GeneralContactFvPatchVectorField::solid4GeneralContactFvPatchVectorField
(
    const solid4GeneralContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
	//******** based on solid General***************
    firstPatchPtr_(NULL),
    firstPatchIndexInListPtr_(NULL),
	currentMasterPtr_(NULL),
    currentMasterIndexPtr_(NULL),
	currentSlavePtr_(NULL),
	localSlavePtr_(NULL),
	//******** END based on solid General*************
    dict_(ptf.dict_),
    //master_(ptf.master_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    writePointDistanceFields_(ptf.writePointDistanceFields_),
    slavePatchNames_(ptf.slavePatchNames_),
    slavePatchIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    slaveZones_(),
    zoneToZones_(),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    contact_(ptf.contact_),
    contactPerSlave_(),
    scaleFaceTractionsNearDownstreamPatch_
    (
        ptf.scaleFaceTractionsNearDownstreamPatch_
    ),
    scaleTractionFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_),
	//******** based on solid General*************
	curPatchTractionPtr_(NULL),
	bbOffset_(ptf.bbOffset_),
	slaveBbOffsetList_(ptf.bbOffset_)
	//******** END based on solid General*************
{
    Info<<"Does it enter here? in C5(ptf, iF)"<<__LINE__<<endl;
	// Do not copy pointer objects
	
	//******************* based on solid General******************
	/*// Copy pointer objects
	if (ptf.firstPatchPtr_)
    {
        firstPatchPtr_ = new bool(*ptf.firstPatchPtr_);
    }

    if (ptf.firstPatchIndexInListPtr_)
    {
        firstPatchIndexInListPtr_ = new label(*ptf.firstPatchIndexInListPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }
	
	if (ptf.curPatchTractionPtr_)
    {
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }
	*/
	//******************* END based on solid General******************
}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::solid4GeneralContactFvPatchVectorField::
~solid4GeneralContactFvPatchVectorField()
{
	Info<<"Does it enter destructor?"<<__LINE__<<endl;
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solid4GeneralContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
	Info<<"Does it call the autoMap(m) in line:"<<__LINE__<<endl;
    if (debug)
    {
        InfoIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::autoMap\n"
            "(\n"
            "    const fvPatchFieldMapper& m\n"
            ")"
        )   << nl << "autoMap: field = " << dimensionedInternalField().name()
            << ", patch = " << patch().name() << endl;
    }

    solidTractionFvPatchVectorField::autoMap(m);

    contact_.autoMap(m);
    scaleTractionFieldPtr_.clear();

    if (contactPerSlave_.size())
    {
        forAll(contactPerSlave_, shadI)
        {
            contactPerSlave_[shadI].autoMap(m);
        }
    }
	
    if (slavePatchNames_)  // if (slavePatchNames_.size() > 0)
    {
        // Let the contact models know about the mapping
        // Be careful, we must pass slave
        // FIX PC 21-Sep-17: move this check inside if (slavePatchNames ... )
        // REPLACE This Check If needed
		if (!firstPatchInList())  // if (!localSlave()) 
        {
            normalModelForThisSlave().autoMap(m);
            frictionModelForThisSlave().autoMap(m);
        }
    }

    // Force all data to be re-created when needed
    clearOut();

    // Reset normal pelanty factors to reinitialise normal models
    normalPenaltyFactors_ =  -1;
}


void Foam::solid4GeneralContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
	Info<<"Does it call the rmap(ptf, addr) in line:"<<__LINE__<<endl;
    solidTractionFvPatchVectorField::rmap(ptf, addr);

    const solid4GeneralContactFvPatchVectorField& dmptf =
        refCast<const solid4GeneralContactFvPatchVectorField>(ptf);

    // PC, I'm not sure if this "if" check should be here...
    if (slavePatchNames_)  //if (slavePatchNames_.size() > 0)
    {
        contact_.rmap(dmptf.contact_, addr);

        if (dmptf.contactPerSlave_.size())
        {
            // Force contactPerSlave to be initialised
            contactPerSlave();

            forAll(contactPerSlave_, shadI)
            {
                contactPerSlave_[shadI].rmap
                (
                    dmptf.contactPerSlave_[shadI], addr
                );
            }
        }
    }

    scaleTractionFieldPtr_.clear();
}

//*************** ORG = NOT based on solidGeneral*****************
//const Foam::wordList&
//*************** END based on solidGeneral*****************
const Foam::List<Foam::word>&
Foam::solid4GeneralContactFvPatchVectorField::slavePatchNames() const
{
	//Info<<"In slavePatchNames() line:"<<__LINE__<<endl;
    if (!slavePatchNames_) //if (slavePatchNames_.size() == 0)
    {
		//Info<<"In slavePatchNames() line:"<<__LINE__<<endl;
       makeSlavePatchNames(); //makeSlavePatchNames(dict_);
    }
	//Info<<"*slavePatchNames_ in slavePatchNames(): "<<*slavePatchNames_<<endl;

    return *slavePatchNames_;
}


const Foam::labelList&
Foam::solid4GeneralContactFvPatchVectorField::slavePatchIndices() const
{
    if (!slavePatchIndicesPtr_)
    {
        makeSlavePatchNames(); //calcSlavePatchIndices();
    }
	//Info<<"*slavePatchIndicesPtr_ in slavePatchIndices(): "<<*slavePatchIndicesPtr_<<endl;

    return *slavePatchIndicesPtr_;
}


const Foam::solid4GeneralContactFvPatchVectorField&
Foam::solid4GeneralContactFvPatchVectorField::slavePatchField() const
{
	Info<<"In slavePatchField() line:"<<__LINE__<<endl;
    /*
	if (slavePatchIndices().size() != 1)
    {
        FatalErrorIn
        (
            "const Foam::solid4GeneralContactFvPatchVectorField&\n"
            "Foam::solid4GeneralContactFvPatchVectorField::slavePatchField() const"
        )   << "This function can only be called for a patch with 1 slave "
            << "patch; this patch has " << slavePatchIndices().size()
            << " slave patches!" << abort(FatalError);
    }
	*/

    const volVectorField& field =
        db().lookupObject<volVectorField>(dimensionedInternalField().name());

    return
        refCast<const solid4GeneralContactFvPatchVectorField>
        (
            field.boundaryField()[slavePatchIndices()[0]]
        );
}


//*******************Start a shadI dependent function *******************
//Foam::PtrList<Foam::generalNormalContactModel>&
Foam::generalNormalContactModel&
Foam::solid4GeneralContactFvPatchVectorField::normalContactModel(const label shadowI)
{
	if (firstPatchInList())
    {
		#if(!normalModelDEBUG)
		Info<<"In normalContactModels() line:"<<__LINE__<<endl;
        #endif
		
		if (normalModels_.size() == 0)
        {
			Info<<"In normalModels() line:"<<__LINE__<<endl;
            makeNormalModels(dict_);
        }

        return normalModels_[shadowI];
    }
    else
    {
		#if(!normalModelDEBUG)
		Info<<"ELSE in normalContactModels() line:"<<__LINE__<<endl;
        #endif
		
		const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );
		
		Info<<"slavePatchIndices()[0] in normalModels() :"<<slavePatchIndices()[0]<<endl;
        solid4GeneralContactFvPatchVectorField& slavePatchField =
            const_cast<solid4GeneralContactFvPatchVectorField&>
            (
                refCast<const solid4GeneralContactFvPatchVectorField>
                (
                    field.boundaryField()[slavePatchIndices()[0]]
                )
            );

        return slavePatchField.normalContactModel(shadowI);
    }
	
	//**********************END
	
	
	/*
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("normalModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (normalModels_.empty())
    {
        calcNormalContactModels(dict_);
    }

    return normalModels_[shadowI];
	*/
}

//const Foam::PtrList<Foam::generalNormalContactModel>&
const Foam::generalNormalContactModel&
Foam::solid4GeneralContactFvPatchVectorField::normalContactModel
(
    const label shadowI
) const
{
	
	if (firstPatchInList())
    {
        if (normalModels_.size() == 0)
        {
            makeNormalModels(dict_);
        }

        return normalModels_[shadowI];
    }
    else
    {
		#if(normalModelDEBUG)
		Info<<"ELSE in normalModels() line:"<<__LINE__<<endl;
        #endif
		
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

        return slavePatchField.normalContactModel(shadowI);
    }
	
	/*
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("normalModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (normalModels_.empty())
    {
        calcNormalContactModels(dict_);
    }

    return normalModels_[shadowI];
	*/
}

//*******************End shadI dependent function **********************


Foam::PtrList<Foam::generalNormalContactModel>&
Foam::solid4GeneralContactFvPatchVectorField::normalModels()
{
    if (firstPatchInList())
    {
		#if(!normalModelDEBUG)
		Info<<"In normalModels() line:"<<__LINE__<<endl;
        #endif
		
		if (normalModels_.size() == 0)
        {
			Info<<"In normalModels() line:"<<__LINE__<<endl;
            makeNormalModels(dict_);
        }

        return normalModels_;
    }
    else
    {
		#if(!normalModelDEBUG)
		Info<<"ELSE in normalModels() line:"<<__LINE__<<endl;
        #endif
		
		const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );
		
		Info<<"slavePatchIndices()[0] in normalModels() :"<<slavePatchIndices()[0]<<endl;
        solid4GeneralContactFvPatchVectorField& slavePatchField =
            const_cast<solid4GeneralContactFvPatchVectorField&>
            (
                refCast<const solid4GeneralContactFvPatchVectorField>
                (
                    field.boundaryField()[slavePatchIndices()[0]]
                )
            );

        return slavePatchField.normalModels();
    }
}


const Foam::PtrList<Foam::generalNormalContactModel>&
Foam::solid4GeneralContactFvPatchVectorField::normalModels() const
{
    if (firstPatchInList())
    {
        if (normalModels_.size() == 0)
        {
            makeNormalModels(dict_);
        }

        return normalModels_;
    }
    else
    {
		#if(normalModelDEBUG)
		Info<<"ELSE in normalModels() line:"<<__LINE__<<endl;
        #endif
		
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

        return slavePatchField.normalModels();
    }
}


Foam::PtrList<Foam::generalFrictionContactModel>&
Foam::solid4GeneralContactFvPatchVectorField::frictionModels()
{
    if (firstPatchInList())
    {
        if (frictionModels_.size() == 0)
        {
            makeFrictionModels(dict_);
        }

        return frictionModels_;
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

        return slavePatchField.frictionModels();
    }
}


const Foam::PtrList<Foam::generalFrictionContactModel>&
Foam::solid4GeneralContactFvPatchVectorField::frictionModels() const
{
    if (firstPatchInList())
    {
        if (frictionModels_.size() == 0)
        {
            makeFrictionModels(dict_);
        }

        return frictionModels_;
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

        return slavePatchField.frictionModels();
    }
}


Foam::generalNormalContactModel&
Foam::solid4GeneralContactFvPatchVectorField::normalModelForThisSlave()
{	
	Info<<"In normalModelForThisSlave() line: "<<__LINE__<<endl;
	#if(normalModelDEBUG)
	Info<< "patch().name() in normalModelForThisSlave() "<<patch().name()<<endl;
	Info<< "patch().index() in normalModelForThisSlave() "<<patch().index()<<endl;
    Info<< "*localSlavePtr_ in normalModelForThisSlave() "<<*localSlavePtr_<<endl;
	#endif
	if (firstPatchInList())
    {
        FatalErrorIn
        (
            "generalNormalContactModel&"
            "solid4GeneralContactFvPatchVectorField::normalModelForThisSlave()"
        )   << "The master is not allowed to called this fucntion!"
            << abort(FatalError);
    }

    // Lookup the master patch corresponding to the current slave patch
    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    if (returnReduce(slavePatchIndices().size() == 0, maxOp<bool>()))
    {
        FatalError
            << "slavePatchIndices().size() == 0" << exit(FatalError);
    }
	
	Info<<"slavePatchIndices()[0] in normalModelForThisSlave(): "<<slavePatchIndices()[0]<<endl;
    Info<<"In normalModelForThisSlave() line: "<<__LINE__<<endl;	
	
	const solid4GeneralContactFvPatchVectorField& masterPatchField =
        refCast<const solid4GeneralContactFvPatchVectorField>
        (
            field.boundaryField()[slavePatchIndices()[0]]
        );
		
	Info<<"patch().name() in normalModelForThisSlave(): "<<patch().name()<<endl;
	
    // The master may have multiple slaves so we need to find which model
    // corresponds to the current slave patch
    const wordList& shadPatchNames = masterPatchField.slavePatchNames();
    label masterSlaveID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterSlaveID = shadPatchI;
			Info<<"In normalModelForThisSlave() line: "<<__LINE__<<endl;
			Info<<"masterSlaveID in normalModelForThisSlave(): "<<masterSlaveID<<endl;
            break;
        }
    }
	
	Info<<"In normalModelForThisSlave() line: "<<__LINE__<<endl;
	Info<< "patch().name() in normalModelForThisSlave() "<<patch().name()<<endl;
	Info<< "patch().index() in normalModelForThisSlave() "<<patch().index()<<endl;
	

    if (masterSlaveID == -1)
    {
        FatalErrorIn
        (
            "void solid4GeneralContactFvPatchVectorField::"
            "normalModelForThisSlave()"
        )   << "Something went wrong when looking for the slavePatch"
            << abort(FatalError);
    }

    return normalModels()[masterSlaveID];
}


Foam::generalFrictionContactModel&
Foam::solid4GeneralContactFvPatchVectorField::frictionModelForThisSlave()
{
    if (firstPatchInList())
    {
        FatalErrorIn
        (
            "generalFrictionContactModel& "
            "solid4GeneralContactFvPatchVectorField::frictionModelForThisSlave()"
        )   << "The master is not allowed to called this function!"
            << abort(FatalError);
    }

    // Lookup the master patch corresponding to the current slave patch
    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solid4GeneralContactFvPatchVectorField& masterPatchField =
        refCast<const solid4GeneralContactFvPatchVectorField>
        (
            field.boundaryField()[slavePatchIndices()[0]]
        );

    // The master may have multiple slaves so we need to find which model
    // corresponds to the current slave patch
    const wordList& shadPatchNames = masterPatchField.slavePatchNames();
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
            "void solid4GeneralContactFvPatchVectorField::"
            "frictionModelForThisSlave()"
        )   << "Something went wrong when looking for the slavePatch"
            << abort(FatalError);
    }

    return frictionModels()[masterSlaveID];
}


void Foam::solid4GeneralContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
		Info<<"if (updated()) check in solid4GeneralContact's updateCoeffs()"<<__LINE__<<endl;
        return;
    }
	#if(!ISDEBUG)
	Info<< "Check 1: patch().name() in updateCoeffs() "<<patch().name()<<endl;
	Info<< "patch().index() in updateCoeffs() "<<patch().index()<<endl;	
	Info<< "slavePatchNames().size() in updateCoeffs() "<<slavePatchNames().size()<<endl;	
	
	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;	
	#endif
	//*************** based on solidGeneral*****************
	boolList activeContactPairs(slavePatchNames().size(), false);  // false);
	//*************** END based on solidGeneral**************
	#if(ISDEBUG)
	Info<<"activeContactPairs in updateCoeffs() "<<activeContactPairs<<endl;
	
	Info<< "this->db().time().timeIndex() "<<this->db().time().timeIndex()<<endl;
	Info<< "curTimeIndex_ "<<curTimeIndex_<<endl;
    #endif
	
	// Only the local masters calculates the contact force and the local
        // master interpolates this force
        const boolList& locSlave = localSlave();
		
	if (curTimeIndex_ != this->db().time().timeIndex())
    {
		Info<<"In updateCoeffs() line:"<<__LINE__<<endl;		
        // Update old quantities at the start of a new time-step
        curTimeIndex_ = this->db().time().timeIndex();
		

        //if (firstPatchInList())
        //{
			
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
            // Let the contact models know that it is a new time-step, in case
            // they need to update anything
			#if(ISDEBUG)
			Info<<"Is it activeContactPairs? in updateCoeffs() "<<activeContactPairs<<endl;
            #endif
			forAll(activeContactPairs, shadPatchI)
            {
				
				if (locSlave[shadPatchI])
				{	
				Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
                normalModels()[shadPatchI].newTimeStep();
                frictionModels()[shadPatchI].newTimeStep();

                // Force N^2 contact search at least once per time-step
                zoneToZones()[shadPatchI].clearPrevCandidateMasterNeighbors();
				}
			}
        //}
		
    }
	#if(ISDEBUG)
	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
    #endif
	// Move the master and slave zone to the deformed configuration
    if (currentMaster())   //(firstPatchInList())
    {
	#if(ISDEBUG)
	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
	#endif
	moveZonesToDeformedConfiguration();
	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
	}
	#if(ISDEBUG)
	Info<<"zone().patch().localPoints() in updateCoeffs():"<<zone().patch().localPoints()<<endl;
	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
	#endif
	
		
	#if(ISDEBUG)
	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
    #endif
	// Delete the zone-to-zone interpolator weights as the zones have moved
    // const wordList& shadPatchNames = slavePatchNames();
    forAll(activeContactPairs, shadPatchI)
    {
		if (locSlave[shadPatchI])//if (localSlave()[shadPatchI])
		{
		Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
        zoneToZones()[shadPatchI].movePoints
        (
            tensorField(0), tensorField(0), vectorField(0)
        );
		}
    }
	
	// Create master bounding box used for quick check
        boundBox masterBb(zone().patch().localPoints(), false);
		#if(ISDEBUG)
		Info<<"zone().patch().localPoints() in updateCoeffs(): "<<zone().patch().localPoints()<<endl;
		Info<< "patch().name() in updateCoeffs() "<<patch().name()<<endl;
		Info<< "zone().patch().name() in updateCoeffs() "<<zone().patch().name()<<endl;
		Info<<"In updateCoeffs() line:"<<__LINE__<<endl;		
		Info<<"masterBb in updateCoeffs(): "<<masterBb<<endl;
		#endif
		
		// The BB may have zero thickness in one of the directions e.g. for a
        // flat patch, so we will check for this and, if found, create an offset
        const scalar bbOff = bbOffset();
		#if(ISDEBUG)
		Info<<"masterBb.minDim() in updateCoeffs(): "<<masterBb.minDim()<<endl;
		Info<<"What is bbOff? in updateCoeffs():"<<bbOff<<endl;
		#endif
		if (masterBb.minDim() < bbOff)
        {
            const vector bbDiag = masterBb.max() - masterBb.min();
			
			#if(ISDEBUG)
			Info<<"bbDiag in updateCoeffs():"<<bbDiag<<endl;
            #endif
			
			if (bbDiag.x() < bbOff)
            {
				#if(ISDEBUG)
				Info<<"Does it enter this check? in updateCoeffs():"<<__LINE__<<endl;
				Info<<"masterBb.min() in updateCoeffs():"<<masterBb.min()<<endl;
				Info<<"masterBb.max() in updateCoeffs():"<<masterBb.max()<<endl;
                #endif
				
				vector offset(bbOff, 0, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
				#if(ISDEBUG)
				Info<<"What is offset? in updateCoeffs():"<<offset<<endl;				
				#endif
			}
			//else if (bbDiag.y() < bbOff)
            if (bbDiag.y() < bbOff)
            {
				#if(ISDEBUG)
				Info<<"Does it enter this check? in updateCoeffs():"<<__LINE__<<endl;
                #endif
				
				vector offset(0, bbOff, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
				Info<<"What is Y offset? in updateCoeffs():"<<offset<<endl;
            }
			//else if (bbDiag.z() < bbOff)
            if (bbDiag.z() < bbOff)
            {
                vector offset(0, 0, bbOff);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
			#if(ISDEBUG)
			Info<<"masterBb.min() in updateCoeffs():"<<masterBb.min()<<endl;
			Info<<"masterBb.max() in updateCoeffs():"<<masterBb.max()<<endl;
			#endif
		}
		
	// Accumulated traction for the current patch
        vectorField curPatchTraction(patch().size(), vector::zero);
	//	Info<<"What is curPatchTraction? "<<curPatchTraction<<endl;
	
	// The BB may have zero thickness in one of the directions e.g. for a
	// flat patch, so we will check for this and, if found, create an offset
	//	const scalarField slaveBbOff = slaveBbOffset();
	//	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
	//	Info<<"slaveBbOff in updateCoeffs() :"<<slaveBbOff<<endl;
			
	forAll(activeContactPairs,shadPatchI)	
	{
		#if(ISDEBUG)
		Info<< "patch().name() in updateCoeffs() "<<patch().name()<<endl;
		Info<< "patch().index() in updateCoeffs() "<<patch().index()<<endl;
		Info<< "slavePatchNames()[shadPatchI] in updateCoeffs() "<<slavePatchNames()[shadPatchI]<<endl;
		#endif
		Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
		//Info<<"activeContactPairs[shadPatchI] in updateCoeffs():"<<activeContactPairs[shadPatchI]<<endl;
			// Create slave bounding box
            boundBox slaveBb(slaveZones()[shadPatchI].patch().localPoints(), false);
			#if(ISDEBUG)
			Info<<"slaveBb in updateCoeffs(): "<<slaveBb<<endl;
			Info<<"slaveBb.minDim() in updateCoeffs(): "<<slaveBb.minDim()<<endl;
			#endif
			
			Info<<"bbOff in updateCoeffs(): "<<bbOff<<endl;
			// Check for a zero dimension in the slaveBb
            if (slaveBb.minDim() < bbOff)
            {
				#if(ISDEBUG)
				Info<<"In updateCoeffs():"<<__LINE__<<endl;
                #endif
				
				const vector bbDiag = slaveBb.max() - slaveBb.min();

                if (bbDiag.x() < bbOff)
                {
					#if(ISDEBUG)
				Info<<"Does it enter this check? in updateCoeffs():"<<__LINE__<<endl;
				Info<<"slaveBb.min() in updateCoeffs():"<<slaveBb.min()<<endl;
				Info<<"slaveBb.max() in updateCoeffs():"<<slaveBb.max()<<endl;
                    #endif
					
					vector offset(bbOff, 0, 0);
                    slaveBb.min() -= offset;
                    slaveBb.max() += offset;
                }
				//else if (bbDiag.y() < bbOff)
                if (bbDiag.y() < bbOff)
                {
					Info<<"In updateCoeffs():"<<__LINE__<<endl;
                    vector offset(0, bbOff, 0);
                    slaveBb.min() -= offset;
                    slaveBb.max() += offset;
                }
				//else if (bbDiag.z() < bbOff)
                if (bbDiag.z() < bbOff)
                {
					#if(ISDEBUG)
					Info<<"In updateCoeffs():"<<__LINE__<<endl;
                    #endif
					
					vector offset(0, 0, bbOff);
                    slaveBb.min() -= offset;
                    slaveBb.max() += offset;
                }
				
				#if(ISDEBUG)
				Info<<"slaveBb.min() in updateCoeffs():"<<slaveBb.min()<<endl;
				Info<<"slaveBb.max() in updateCoeffs():"<<slaveBb.max()<<endl;
				Info<<"bbDiag in updateCoeffs():"<<bbDiag<<endl;
				Info<<"bbDiag.x() updateCoeffs():"<<bbDiag.x()<<endl;
				Info<<"bbDiag.y() updateCoeffs():"<<bbDiag.y()<<endl;
				Info<<"bbDiag.z() updateCoeffs():"<<bbDiag.z()<<endl;
				#endif
				
			}
			
			#if(ISDEBUG)
			Info<< "patch().name() in updateCoeffs() "<<patch().name()<<endl;
			Info<< "slavePatchNames()[shadPatchI] in updateCoeffs() "<<slavePatchNames()[shadPatchI]<<endl;			
			Info<<"In updateCoeffs():"<<__LINE__<<endl;
			Info<<"masterBb in updateCoeffs():"<<masterBb<<endl;
			Info<<"slaveBb in updateCoeffs():"<<slaveBb<<endl;
			Info<<"masterBb.overlaps(slaveBb) in updateCoeffs():"<<masterBb.overlaps(slaveBb)<<endl;
			#endif
			
			if (masterBb.overlaps(slaveBb))
            {
				Info<<"In updateCoeffs():"<<__LINE__<<endl;
                activeContactPairs[shadPatchI] = true;
            }
			
			#if(ISDEBUG)
			Info<<"activeContactPairs[shadPatchI] in updateCoeffs():"<<activeContactPairs[shadPatchI]<<endl;
			#endif
			
		if (activeContactPairs[shadPatchI])
        {
			if(currentMaster())
			{
				#if(ActivePairDEBUG)
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;	
				#endif
			}
			
			Info<<" "<<endl;
			#if(!ActivePairDEBUG)
			Info<< "patch INDEX in updateCoeffs() "<<patch().index()<<endl;
			Info<<"Checking mMASTER or sSLAVE in updateCoeffs() line:"<<__LINE__<<endl;
			Info<<"locSlave[shadPatchI] in updateCoeffs(): "<<locSlave[shadPatchI]<<endl;
			Info<< "patch().name() in updateCoeffs() "<<patch().name()<<endl;
			Info<<"slavePatchNames in updateCoeffs() "<<slavePatchNames()[shadPatchI]<<endl;
			#endif
			Info<<"...................................... "<<endl;
			
			if (locSlave[shadPatchI])  //MASTER starts
            {
				Info<<"MASTER in updateCoeffs() line:"<<__LINE__<<endl;
			// Reset the traction to zero as we will accumulate it over all the
			// slave patches
			//traction() = vector::zero;
			
			#if(masterOfPairDEBUG)
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
            #endif
			
			// Calculate the slave patch face unit normals as they are used by
            // both the normal and friction models
            const vectorField slavePatchFaceNormals =
                slaveZones()[shadPatchI].globalFaceToPatch
                (
                    slaveZones()[shadPatchI].globalPatch().faceNormals()
                );
			
			#if(masterOfPairDEBUG)
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
            #endif
			
			// Interpolate the master displacement increment to the slave patch
            // as it is required by specific normal and friction contact models

            vectorField patchDD(patch().size(), vector::zero);
            vectorField slavePatchDD
            (
                patch().boundaryMesh()[slavePatchIndices()[shadPatchI]].size(),
                vector::zero
            );
			
			#if(masterOfPairDEBUG)
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
			#endif
			
				if (movingMesh())
				{
                // Updated Lagrangian, we will directly lookup the displacement
                // increment

                const volVectorField& DD =
                    db().lookupObject<volVectorField>("DU");
			//		db().lookupObject<volVectorField>("DD");

                patchDD = DD.boundaryField()[patch().index()];
                slavePatchDD =
                    DD.boundaryField()[slavePatchIndices()[shadPatchI]];
				}
				else
				{
					Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
                // We will lookup the total displacement and old total
                // displacement

                const volVectorField& D =
                    db().lookupObject<volVectorField>("U");
				//	db().lookupObject<volVectorField>("D");
				
				#if(masterOfPairDEBUG)
				Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
				Info<<"patch().index() in updateCoeffs(): "<<patch().index()<<endl;
                #endif
				
				patchDD =
                    D.boundaryField()[patch().index()]
                  - D.oldTime().boundaryField()[patch().index()];
                slavePatchDD =
                    D.boundaryField()[slavePatchIndices()[shadPatchI]]
                  - D.oldTime().boundaryField()
                    [
                        slavePatchIndices()[shadPatchI]
                    ];
				}
			
            // Master zone DD
            const vectorField zoneDD = zone().patchFaceToGlobal(patchDD);
			
			#if(masterOfPairDEBUG)
			Info<< "zoneDD in updateCoeffs() "<<zoneDD<< endl;
			
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
            #endif
			
			// Master patch DD interpolated to the slave patch
            const vectorField patchDDInterpToSlavePatch =
                slaveZones()[shadPatchI].globalFaceToPatch
                (
                    zoneToZones()[shadPatchI].masterToSlave(zoneDD)()
                );
			
			#if(masterOfPairDEBUG)
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
            #endif
			
			// Calculate normal contact forces
            // slavePatchDD is the DU on the slave patch, whereas
            // patchDDInterpToSlavePatch is the master patch DU interpolated to
            // the slave; and the difference between these two is the slip (and
            // also the normal component of DU)
            normalModels()[shadPatchI].correct
            (
                slavePatchFaceNormals,
                slaveZones()[shadPatchI].globalPointToPatch
                (
                    zoneToZones()[shadPatchI].slavePointDistanceToIntersection()
                ),
                // zoneToZones()[shadPatchI],
                slavePatchDD,
                patchDDInterpToSlavePatch
            );

            // Calculate friction contact forces
            frictionModels()[shadPatchI].correct
            (
                normalModels()[shadPatchI].slavePressure(),
                slavePatchFaceNormals,
                normalModels()[shadPatchI].areaInContact(),
                slavePatchDD,
                patchDDInterpToSlavePatch
            );

				if (rigidMaster_)
				{
                // Set to master to traction free to mimic a rigid contact
                traction() = vector::zero;

                // Set contact indicator field
                contactPerSlave()[shadPatchI] = 0.0;
				}
				else
				{
				// Interpolate slave traction to the master
                const vectorField slavePatchTraction =
                   - frictionModels()[shadPatchI].slaveTractionForMaster()
                   - normalModels()[shadPatchI].slavePressure();

                const vectorField slaveZoneTraction =
                    slaveZones()[shadPatchI].patchFaceToGlobal
                    (
                        slavePatchTraction
                    );

                // We have two options for interpolating from the slave to the
                // master:
                // 1. face-to-face
                // 2. point-to-point
                // We will use 1.

                // Calculate traction for this contact
					vectorField tractionForThisSlave =
                    zone().globalFaceToPatch
                    (
                        zoneToZones()[shadPatchI].slaveToMaster
                        (
                            slaveZoneTraction
                        )()
                    );

                // Accumulate the traction on the master patch
                curPatchTractions(shadPatchI) = tractionForThisSlave;
				//traction() += tractionForThisSlave;
				curPatchTraction += curPatchTractions(shadPatchI);
				
				#if(!masterOfPairDEBUG)
				Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
				Info<<"MASTER curPatchTraction in updateCoeffs() :"<<curPatchTraction<<endl;
                #endif
				
				// Update contactPerSlave field
                // Note: this is used by thermalContact to know which faces
                // are in contact
                const scalarField magTraction = mag(tractionForThisSlave);
                const scalar tol = 1e-6*gMax(magTraction);
                scalarField& contactForThisSlave =
                    contactPerSlave()[shadPatchI];
					forAll(contactForThisSlave, faceI)
					{
						if (magTraction[faceI] > tol)
						{
                        contactForThisSlave[faceI] = 1.0;
						}
						else
						{
                        contactForThisSlave[faceI] = 0.0;
						}
					}
				}
				
				#if(!masterOfPairDEBUG)
			Info<<"End of MASTER computation in updateCoeffs() line:"<<__LINE__<<endl;
			Info<<"Which pair of this - MASTER in updateCoeffs(): "<<shadPatchI<<endl;
				#endif
			Info<<"....................end........... "<<endl;	
			}
			else
			{
			
			Info<<"SLAVE in updateCoeffs() line:"<<__LINE__<<endl;
			
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
                                slavePatchIndices()[shadPatchI]
                            ]
                        );
						
			
			const label masterShadowI =
                        localMasterField.findSlaveID(patch().index());
        
		//****************** Start Test with shadPatchI dependent function **************
		
		curSlaveTractions(shadPatchI)= localMasterField.normalContactModel(masterShadowI).slavePressure();
			
		
		//****************** End test with shadPatchI dependent function **************
		
		// Set the traction on the slave patch
        // The master stores the friction and normal models, so we need to find
        // which models correspond to the current slave
        //traction() =
		curPatchTractions(shadPatchI) =
            frictionModelForThisSlave().slaveTraction()
          + normalModelForThisSlave().slavePressure();
		
		Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
		
		Info<<"shadPatchI in updateCoeffs(): "<<shadPatchI<<endl;
		Info<<"slavePatchIndices()[shadPatchI] in updateCoeffs(): "<<slavePatchIndices()[shadPatchI]<<endl;
		
		curPatchTraction += curPatchTractions(shadPatchI);
		
		Info<<"SLAVE curPatchTraction in updateCoeffs() :"<<curPatchTraction<<endl;
		
        // TESTING - START
        // Scale traction vectors on faces, which share an edge with the
        // downstream patch
        // This is an attempt to fix an issue where the first row of faces
        // deform unphysically when being drawn into the die
				if (scaleFaceTractionsNearDownstreamPatch_)
				{
					Info<<"TESTING - START in updateCoeffs() line:"<<__LINE__<<endl;
				//traction() *= scaleTractionField();
				curPatchTractions(shadPatchI) *= scaleTractionField();
				}
        // TESTING - END
		Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
		
        // Update contactPerSlave field
        // Note: this is used by thermalContact to know which faces
        // are in contact
        //const scalarField magTraction = mag(traction());
        const scalarField magTraction = mag(curPatchTractions(shadPatchI));
		const scalar tol = 1e-6*gMax(magTraction);
        scalarField& contactForThisSlave = contactPerSlave()[0];
				forAll(contactForThisSlave, faceI)
				{
					Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
					if (magTraction[faceI] > tol)
					{
					contactForThisSlave[faceI] = 1.0;
					}
					else
					{
					contactForThisSlave[faceI] = 0.0;
					}
				}
				
				Info<<"slavePatchNames()[shadPatchI] in updateCoeffs(): "<<slavePatchNames()[shadPatchI]<<endl;
				Info<<"slavePatchIndices()[shadPatchI] in updateCoeffs(): "<<slavePatchIndices()[shadPatchI]<<endl;
				Info<<"Which pair of this - SLAVE in updateCoeffs(): "<<shadPatchI<<endl;
				Info<<"End of SLAVE computation in updateCoeffs(): "<<__LINE__<<endl;
				Info<<".................end.............. "<<endl;				
			}// SLAVE \n
		} // if contact pair is active
	}// forAll contact pairs
	
	Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
	
	// Set master gradient based on accumulated traction
    traction() = curPatchTraction;
	
	Info<<"traction() in updateCoeffs():"<<traction()<<endl;

    // Accumulate the contact indicator field
    contact_ = 0.0;
    PtrList<scalarField>& contactPerSlave = this->contactPerSlave();
		forAll(contactPerSlave, shadI)
		{
		Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
        contact_ += contactPerSlave[shadI];
		}

    // Scale any face in contact with more than one slave
	//******************* START Scaling any face in contact ******************
     	if (gMax(contact_) > (1.0 + SMALL))
		{
		Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
			forAll(contact_, faceI)
			{
				if (contact_[faceI] > (1.0 + SMALL))
				{
				Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
                // Update the contact weights corresponding to each slave
                scalar sumContact = 0.0;
					forAll(contactPerSlave, shadI)
					{
                    contactPerSlave[shadI][faceI] /= contact_[faceI];
                    sumContact += contactPerSlave[shadI][faceI];
					}

					if (sumContact > (1.0 + SMALL))
					{
                    FatalErrorIn
                    (
                        "void solid4GeneralContactFvPatchVectorField::"
                        "updateCoeffs()"
                    )   << "There is a problem normalising the contact field"
                        << ", sumContact is: " << sumContact
                        << abort(FatalError);
					}

                // Reset accumulated contact face value to 1.0
                contact_[faceI] = 1.0;
				}
			}
		} 
	//******************* END Scaling any face in contact ******************
	Info<<"Before solidTractionFvPatch in updateCoeffs() line:"<<__LINE__<<endl;
    solidTractionFvPatchVectorField::updateCoeffs();
}

/*
Foam::tmp<Foam::scalarField>
Foam::solid4GeneralContactFvPatchVectorField::frictionHeatRate() const
{
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
    // Consider storing frictionHeatRate instead of recalculating multiple times

    if (!currentMaster())//(!firstPatchInList())
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::scalarField> Foam::"
            "solid4GeneralContactFvPatchVectorField::frictionHeatRate() const"
        )   << "Only master can call frictionHeatRate function!"
            << abort(FatalError);
    }
	
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

    // For now, we assume traction is constant over time-step
    // Todo: we should use trapezoidal rule
    vectorField curTraction(patch().size(), vector::zero);
	
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

    tmp<scalarField> tfrictionHeatRate
    (
        new scalarField(curTraction.size(), 0.0)
    );
    scalarField& frictionHeatRate = tfrictionHeatRate();
	
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

    forAll(slavePatchNames(), shadPatchI)
    {
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
        // Calculate slip

        const vectorField slavePatchSlip = frictionModels()[shadPatchI].slip();
		
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
		
        const vectorField slaveZoneSlip =
            slaveZones()[shadPatchI].patchFaceToGlobal
            (
                slavePatchSlip
            );
			
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

        // Interpolate from slave to master

        const vectorField masterZoneSlip =
            zoneToZones()[shadPatchI].slaveToMaster(slaveZoneSlip);
			
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

        const vectorField masterPatchSlip =
            zone().globalFaceToPatch
            (
                masterZoneSlip
            );

        const scalar deltaT =
            patch().boundaryMesh().mesh().time().deltaTValue();

        // Accumulate frictionHeatRate for each slave patch

        // Rate of dissipated frictional energy for this timestep
        // The dot product of the traction vectors and the slip vectors gives
        // the dissipated frictional energy per unit area; which is always
        // positive
        frictionHeatRate += mag(traction() & (masterPatchSlip/deltaT));
    }
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

    return tfrictionHeatRate;
}
*/

// ***************** Based on solidGeneral *******************
Foam::tmp<Foam::scalarField>
Foam::solid4GeneralContactFvPatchVectorField::frictionHeatRate() const
{
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
    // Consider storing frictionHeatRate instead of recalculating multiple times

    if (!currentMaster())//(!firstPatchInList())
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::scalarField> Foam::"
            "solid4GeneralContactFvPatchVectorField::frictionHeatRate() const"
        )   << "Only master can call frictionHeatRate function!"
            << abort(FatalError);
    }
	
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

    // For now, we assume traction is constant over time-step
    // Todo: we should use trapezoidal rule
    vectorField curTraction(patch().size(), vector::zero);
	
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

    tmp<scalarField> tfrictionHeatRate
    (
        new scalarField(curTraction.size(), 0.0)
    );
    scalarField& frictionHeatRate = tfrictionHeatRate();
	
	const scalar deltaT =
            patch().boundaryMesh().mesh().time().deltaTValue();
	
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
	
	const boolList& locSlave = localTSlave();

    forAll(locSlave, shadPatchI) //forAll(slavePatchNames(), shadPatchI)
    {
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
       
		//**************** based on solidGeneral *****************
		vectorField curPatchSlip(patch().size(), vector::zero);
		
		// Calculate slip
		if (locSlave[shadPatchI])
        {
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
			const vectorField slavePatchSlip = frictionModels()[shadPatchI].slip();
			
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
		
		Info<<"patch().index() in frictionHeatRate(): "<<patch().index()<<endl;
		Info<<"patch().name() in frictionHeatRate(): "<<patch().name()<<endl;
		
        const vectorField slaveZoneSlip =
            slaveZones()[shadPatchI].patchFaceToGlobal
            (
                slavePatchSlip
            );
			
			
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

        // Interpolate from slave to master

        const vectorField masterZoneSlip =
            zoneToZones()[shadPatchI].slaveToMaster(slaveZoneSlip);
			
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

        //const vectorField masterPatchSlip =
            curPatchSlip =
			zone().globalFaceToPatch
            (
                masterZoneSlip
            );
		
		/*const scalar deltaT =
            patch().boundaryMesh().mesh().time().deltaTValue();
			*/
		
		}
		else 
		{
		Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;
			curPatchSlip = frictionModels()[shadPatchI].slip();
		//const vectorField slavePatchSlip = frictionModels()[shadPatchI].slip();
		}
		
        // Accumulate frictionHeatRate for each slave patch

        // Rate of dissipated frictional energy for this timestep
        // The dot product of the traction vectors and the slip vectors gives
        // the dissipated frictional energy per unit area; which is always
        // positive
		frictionHeatRate += mag(traction() & (curPatchSlip/deltaT));
        //frictionHeatRate += mag(traction() & (masterPatchSlip/deltaT));
		//**************** End based on solidGeneral *****************
    
	}	
	
	Info<<"In frictionHeatRate() line:"<<__LINE__<<endl;

    return tfrictionHeatRate;
}
// ***************** End Based on solidGeneral *****************

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


void Foam::solid4GeneralContactFvPatchVectorField::write(Ostream& os) const
{
	Info<<"In solid4GeneralContact::write function "<<__LINE__<<endl;
    // If the slavePatchIndices pointer is not set then we will assume that the
    // contact models were not created and nothing has changed; so we will just
    // output the input dict unchanged
    if (!slavePatchNames_)  //if (slavePatchNames_.size() == 0)
    {
		Info<<"In solid4GeneralContact::write function "<<__LINE__<<endl;
        // Overwrite fields in the dict
        dictionary& dict = const_cast<dictionary&>(dict_);

        dict.remove("gradient");
        dict.remove("value");
        dict.remove("traction");
        dict.remove("pressure");

        //dict.add("gradient", gradient());
        const vectorField& patchValue = *this;
        //dict.add("value", patchValue);
        //dict.add("traction", traction());
        //dict.add("pressure", pressure());

        // Write the dictionary
        dict_.write(os, false);

        gradient().writeEntry("gradient", os);
        patchValue.writeEntry("value", os);
        traction().writeEntry("traction", os);
        pressure().writeEntry("pressure", os);

        return;
    }
	
		if(
	    !localSlavePtr_
		&& dimensionedInternalField().name() == "U_0"
		)
		{
//	  	Info<<"Here I am in first U_0 check in write()"<<__LINE__<<endl;
		return;
		} 
	
    solidTractionFvPatchVectorField::write(os);
	Info<<"In solid4GeneralContact::write function "<<__LINE__<<endl;
	/* os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;
	*/	
    const wordList& shadPatchNames = slavePatchNames();
    if (shadPatchNames.size() == 1)
    {		
        os.writeKeyword("slavePatch")
            << shadPatchNames[0] << token::END_STATEMENT << nl;
    }
    else
    {		
        slavePatchNames().writeEntry("slavePatches", os);
	}		
    
	os.writeKeyword("regionOfInterest")
        << regionOfInterest_ << token::END_STATEMENT << nl;
    os.writeKeyword("regionOfInterestTopCorner")
        << regionOfInterestTopCorner_ << token::END_STATEMENT << nl;
    os.writeKeyword("regionOfInterestBottomCorner")
        << regionOfInterestBottomCorner_ << token::END_STATEMENT << nl;
    os.writeKeyword("writeZoneVTK")
        << writeZoneVTK_ << token::END_STATEMENT << nl;
    os.writeKeyword("writePointDistanceFields")
        << writePointDistanceFields_ << token::END_STATEMENT << nl;
    os.writeKeyword("scaleFaceTractionsNearDownstreamPatch")
        << scaleFaceTractionsNearDownstreamPatch_ << token::END_STATEMENT << nl;
    if (scaleFaceTractionsNearDownstreamPatch_)
    {
        os.writeKeyword("downstreamScaleFactor")
            << readScalar(dict_.lookup("downstreamScaleFactor"))
            << token::END_STATEMENT << nl;
        os.writeKeyword("downstreamPatchName")
            << word(dict_.lookup("downstreamPatchName")) << token::END_STATEMENT
            << nl;
    }
	Info<<"In solid4GeneralContact::write function "<<__LINE__<<endl;
    //if (firstPatchInList())   //if (master_) //if (localSlave()[shadPatchI]) 
    //const boolList& locSlave = localSlave();
	//forAll(activeContactPairs, shadPatchI)
	//{
	//if (locSlave[shadPatchI])
    //{
	//if (localSlave()[shadPatchI]) 
	{
        os.writeKeyword("rigidMaster") << rigidMaster_
            << token::END_STATEMENT << nl;
			
		
	// Write the dict from the first contact model

    const label slaveI = 0;
	
	Info<<"In solid4GeneralContact::write function "<<__LINE__<<endl;
	
	Info<< "The current field in solid4GeneralContact::write function is "<< dimensionedInternalField().name()<< endl;
	Info<< "The current patch in solid4GeneralContact::write function is "<< patch().name()<< endl;	
	
	//Info<<"*localSlavePtr_ in solid4GeneralContact::write function "<<*localSlavePtr_<<endl;
	
    if(!localSlavePtr_) //remove this check later, since localSlave should re-compute the local slave
        FatalError  << "solid4GeneralContactFvPatchVectorField::write: localSlavePtr_ NOT defined:" 
                    << "Cannot write slave information because no slave identified!"  
                    << exit(FatalError);; 
					

    if (localSlave()[slaveI])
    {
        os.writeKeyword("generalNormalContactModel")
            << normalModels()[slaveI].type() << token::END_STATEMENT << nl;
        normalModels()[slaveI].writeDict(os);

        os.writeKeyword("generalFrictionContactModel")
            << frictionModels()[slaveI].type() << token::END_STATEMENT << nl;
        frictionModels()[slaveI].writeDict(os);
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
                    slavePatchIndices()[slaveI]
                ]
            );
	/*	
	//My code test
	{
		const labelList slaveIDs = slavePatchIndices();
		label shadPatchI = -1;
		forAll(slaveIDs, I)
		{
			if (patch().index() == slaveIDs[I])
			{
				Info<<"In solid4GeneralContact::write() line:"<<__LINE__<<endl;
				shadPatchI = I;
				break;
			}
			Info<<"In solid4GeneralContact::write() line:"<<__LINE__<<endl;
		}
		Info<<"In solid4GeneralContact::write() line:"<<__LINE__<<endl;

		if (shadPatchI == -1)
		{
			FatalErrorIn("solid4GeneralContact::write(os)")
				<< "slave patch not found!" << abort(FatalError);
		}
	}
	*/
        const label localSlaveID =
            localSlaveField.findSlaveID(patch().index());

        os.writeKeyword("generalNormalContactModel")
            << localSlaveField.normalModels()[localSlaveID].type()
            << token::END_STATEMENT << nl;
        localSlaveField.normalModels()[localSlaveID].writeDict(os);

        os.writeKeyword("generalFrictionContactModel")
            << localSlaveField.frictionModels()[localSlaveID].type()
            << token::END_STATEMENT << nl;
        localSlaveField.frictionModels()[localSlaveID].writeDict(os);
    }
		
        /*
		if (slavePatchNames_.size() == 1)
        {
            os.writeKeyword("generalNormalContactModel")
                << normalModels()[0].type()
                << token::END_STATEMENT << nl;
            normalModels()[0].writeDict(os);

            os.writeKeyword("generalFrictionContactModel")
                << frictionModels()[0].type()
                << token::END_STATEMENT << nl;
            frictionModels()[0].writeDict(os);

            os.writeKeyword("useNewPointDistanceMethod")
                << dict_.lookupOrDefault<Switch>
                (
                    "useNewPointDistanceMethod", false
                )
                << token::END_STATEMENT << nl;

            os.writeKeyword("projectPointsToPatchBoundary")
                << dict_.lookupOrDefault<Switch>
                (
                    "projectPointsToPatchBoundary", false
                )
                << token::END_STATEMENT << nl;

            os.writeKeyword("checkPointDistanceOrientations")
                << dict_.lookupOrDefault<Switch>
                (
                    "checkPointDistanceOrientations", false
                )
                << token::END_STATEMENT << nl;

            os.writeKeyword("usePrevCandidateMasterNeighbors")
                << dict_.lookupOrDefault<Switch>
                (
                    "usePrevCandidateMasterNeighbors", false
                )
                << token::END_STATEMENT << nl;
        }
		*/		
        //else
        //{ 
		/*
		boolList activeContactPairs(slavePatchNames().size(), true);
		const boolList& locSlave = localSlave();
		forAll(activeContactPairs, shadPatchI)
		{
			if (activeContactPairs[shadPatchI])
			{
			Info<<"In updateCoeffs() line:"<<__LINE__<<endl;
				if (locSlave[shadPatchI])
				os  << patch().name() << "_to_"
			
			}
		}	
		*/
			/*
			wordList& slavePatchNames = *slavePatchNames_;
            Info<<"slavePatchNames in solid4GeneralContact::write function "<<slavePatchNames<<endl;
			Info<<"patch().name() in solid4GeneralContact::write function "<<patch().name()<<endl;
			Info<<"In solid4GeneralContact::write function "<<__LINE__<<endl;
			forAll(slavePatchNames, shadPatchI)
            {
				Info<<"In solid4GeneralContact::write function "<<__LINE__<<endl;
                os  << patch().name() << "_to_"
                    << slavePatchNames_[shadPatchI] << "_dict" << nl
                    << '{' << endl;
				Info<<"In solid4GeneralContact::write function "<<__LINE__<<endl;
				
                os.writeKeyword("generalNormalContactModel")
                    << normalModels()[shadPatchI].type()
                    << token::END_STATEMENT << nl;
                normalModels()[shadPatchI].writeDict(os);

                os.writeKeyword("generalFrictionContactModel")
                    << frictionModels()[shadPatchI].type()
                    << token::END_STATEMENT << nl;
                frictionModels()[shadPatchI].writeDict(os);

                os  << '}' << endl;
            }
			*/
        //}
    }
	
    if (writeZoneVTK_)
    {
        if
        (
            dimensionedInternalField().name() == "D"
         || dimensionedInternalField().name() == "DD"
        )
        {
            Info<< "Writing deformed zones to VTK" << endl;
            const word timeName =
                patch().boundaryMesh().mesh().time().timeName();

            zone().globalPatch().writeVTK("zone_" + timeName);

            forAll(slaveZones(), shadI)
            {
                slaveZones()[shadI].globalPatch().writeVTK
                (
                    "slaveZone_" + timeName
                );
            }
        }
    }


    // Write out point distance fields for master and slave
    if (writePointDistanceFields_ && master())
    {
        if (normalModels().size() != 1)
        {
            FatalErrorIn
            (
                "void solid4GeneralContactFvPatchVectorField::"
                "write(Ostream& os) const"
            )   << "The 'writePointDistanceFields' is currently only "
                << "implemented for one-to-one contact"
                << abort(FatalError);
        }

        // Take a reference to the mesh for convenience
        const polyMesh& mesh = patch().patch().boundaryMesh().mesh();

        // Create the point mesh, which is needed for the point field
        pointMesh pMesh(mesh);

        // Create the point distance fields

        pointScalarField dist
        (
            IOobject
            (
                "pointDistance",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        pointVectorField distVecs
        (
            IOobject
            (
                "pointDistanceVectors",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("zero", dimless, vector::zero)
        );

        // Transfer the patch point distances into the dist point field
        {
            // Lookup the master point distance to intersection
            const scalarField masterpd =
                zone().globalPointToPatch
                (
                    zoneToZones()[0].masterPointDistanceToIntersection()
                );
            const vectorField masterpdVecs =
                zone().globalPointToPatch
                (
                    zoneToZones()[0].masterPointDistanceVectorsToIntersection()
                );

            const labelList& masterMeshPoints = patch().patch().meshPoints();

            forAll(masterpd, pI)
            {
                const label pointID = masterMeshPoints[pI];
                dist[pointID] = masterpd[pI];
                distVecs[pointID] = masterpdVecs[pI];
            }
        }

        {
            const scalarField slavepd =
                slaveZones()[0].globalPointToPatch
                (
                    zoneToZones()[0].slavePointDistanceToIntersection()
                );
            const vectorField slavepdVecs =
                slaveZones()[0].globalPointToPatch
                (
                    zoneToZones()[0].slavePointDistanceVectorsToIntersection()
                );

            const labelList& slaveMeshPoints =
                patch().patch().boundaryMesh()
                [
                    slavePatchIndices()[0]
                ].meshPoints();

            forAll(slavepd, pI)
            {
                const label pointID = slaveMeshPoints[pI];
                dist[pointID] = slavepd[pI];
                distVecs[pointID] = slavepdVecs[pI];
            }
        }

        // Write the field
        InfoIn
        (
            "void Foam::solid4GeneralContactFvPatchVectorField::"
            "write(Ostream& os) const"
        )   << "Writing point distance fields: " << dist.name()
            << " and " << distVecs.name() << endl;
        dist.write();
        distVecs.write();
    }
	
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        solid4GeneralContactFvPatchVectorField
    );
}


// ************************************************************************* //
