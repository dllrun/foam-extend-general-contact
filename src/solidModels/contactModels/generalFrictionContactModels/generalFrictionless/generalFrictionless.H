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
    generalFrictionless

Description
    generalFrictionless contact law

SourceFiles
    generalFrictionless.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef generalFrictionless_H
#define generalFrictionless_H

#include "generalFrictionContactModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class generalFrictionless Declaration
\*---------------------------------------------------------------------------*/

class generalFrictionless
:
    public generalFrictionContactModel
{
    // Private data

       //- slave friction traction - always zero for generalFrictionless
       vectorField slaveTraction_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        //generalFrictionless(const generalFrictionless&);

        //- Disallow default bitwise assignment
        void operator=(const generalFrictionless&);

public:

    //- Runtime type information
    TypeName("generalFrictionless");

    // Static data members


    // Constructors

        //- Construct from dictionary
        generalFrictionless
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID
        );

        //- Construct as copy
        generalFrictionless(const generalFrictionless&);

        //- Construct and return a clone
        virtual autoPtr<generalFrictionContactModel> clone() const
        {
            return autoPtr<generalFrictionContactModel>
            (
                new generalFrictionless(*this)
            );
        }


    //- Destructor
    virtual ~generalFrictionless();


    // Member Functions

        //- Correct the friction model
        virtual void correct
        (
            const vectorField& slavePressure,
            const vectorField& slaveFaceNormals,
            const scalarField& slavePatchAreaInContact,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        )
        {}

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return slaveTraction_;
        }

        //- Return slip field
        virtual const vectorField& slip() const
        {
            // slaveTraction is always zero
            return slaveTraction_;
        }

        //- Perform any reauired field mapping e.g. after a topo change
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const
        {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
