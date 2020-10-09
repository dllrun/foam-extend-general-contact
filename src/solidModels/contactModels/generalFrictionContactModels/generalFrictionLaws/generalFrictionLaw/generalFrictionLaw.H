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
    generalFrictionLaw

Description
    Friction law e.g. Coulomb, Orowan, Shaw, etc.
    Given the pressure (or other variable for more complex laws),
    the friction law returns the shear traction.

SourceFiles
    generalFrictionLaw.C
    newFrictionLaw.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef generalFrictionLaw_H
#define generalFrictionLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class generalFrictionContactModel;

/*---------------------------------------------------------------------------*\
                         Class generalFrictionLaw Declaration
\*---------------------------------------------------------------------------*/

class generalFrictionLaw
{
    // Private data

        //- Reference to the generalFrictionContactModel
        const generalFrictionContactModel& frictionModel_;

        //- Friction law dictionary
        dictionary frictionLawDict_;


    // Private Member Functions

        //- Disallow copy construct
        //generalFrictionLaw(const generalFrictionLaw&);

        //- Disallow default bitwise assignment
        void operator=(const generalFrictionLaw&);


public:

    //- Runtime type information
    TypeName("generalFrictionLaw");

    // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            generalFrictionLaw,
            dictionary,
            (
                const word name,
                const generalFrictionContactModel& fricModel,
                const dictionary& dict
            ),
            (name, fricModel, dict)
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<generalFrictionLaw> New
        (
            const word& name,
            const generalFrictionContactModel& fricModel,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        generalFrictionLaw
        (
            const word& name,
            const generalFrictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construct as a copy
        generalFrictionLaw(const generalFrictionLaw& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<generalFrictionLaw> clone() const = 0;


    // Destructor

        virtual ~generalFrictionLaw()
        {}


    // Member Functions

        //- Return reference to friction contact model
        const generalFrictionContactModel& frictionModel() const
        {
            return frictionModel_;
        }

        //- Return friction law dict
        const dictionary& frictionLawDict() const
        {
            return frictionLawDict_;
        }

        //- Return traction that would cause slipping
        virtual scalar slipTraction(const scalar pressure);

        //- Return traction that would cause slipping
        virtual scalar slipTraction
        (
            const scalar pressure, const vector& slipDir
        );

        //- Return traction that would cause slipping
        virtual scalar slipTraction
        (
                const scalar contactPressure,         // Contact pressure
                const vector& faceSlip,               // Slip vector
                const vector& slaveFaceVelocity,      // Velocity of slave face
                const vector& masterFaceVelocity      // Velocity of master face
        );

        //- Return traction that would cause slipping
        virtual scalar slipTraction
        (
            const scalar contactPressure,         // Contact pressure
            const vector& faceSlip,               // Slip vector
            const vector& slaveFaceVelocity,      // Velocity of slave face
            const vector& masterFaceVelocity,     // Velocity of master face
            const label slavePatchIndex,          // Slave patch index
            const label faceIndex                 // Local slave face ID
        );

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
