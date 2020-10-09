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

\*---------------------------------------------------------------------------*/

#include "generalFrictionless.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(generalFrictionless, 0);
    addToRunTimeSelectionTable(generalFrictionContactModel, generalFrictionless, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::generalFrictionless::generalFrictionless
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID
)
:
    generalFrictionContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID
    ),
    slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero)
{}


// Construct as a copy
Foam::generalFrictionless::generalFrictionless(const generalFrictionless& fricLaw)
:
    generalFrictionContactModel(fricLaw),
    slaveTraction_(fricLaw.slaveTraction_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::generalFrictionless::~generalFrictionless()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::generalFrictionless::autoMap(const fvPatchFieldMapper& m)
{
    generalFrictionContactModel::autoMap(m);

    if (debug)
    {
        InfoIn
        (
            "void generalFrictionless::autoMap(const fvPatchFieldMapper& m)"
        )   << "autoMap" << endl;
    }

    slaveTraction_.autoMap(m);
}

// ************************************************************************* //
