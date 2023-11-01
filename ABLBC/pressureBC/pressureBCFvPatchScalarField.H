/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::pressureBCFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a free-stream condition for pressure.

    It is an outlet-inlet condition that uses the velocity orientation to
    continuously blend between zero gradient for normal inlet and fixed value
    for normal outlet flow.

Usage
    \table
        Property        | Description             | Required    | Default value
        U               | velocity field name     | no          | U
        freestreamValue | freestream pressure     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureBC;
        freestreamValue uniform 1e5;
    }
    \endverbatim

Note
    This condition is designed to operate with a freestreamVelocity condition

See also
    Foam::mixedFvPatchField
    Foam::freestreamFvPatchField

SourceFiles
    pressureBCFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureBCFvPatchScalarField_H
#define pressureBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class pressureBCFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class pressureBCFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;


public:

    //- Runtime type information
    TypeName("pressureBC");


    // Constructors

        //- Construct from patch and internal field
        pressureBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureBCFvPatchScalarField onto
        //  a new patch
        pressureBCFvPatchScalarField
        (
            const pressureBCFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureBCFvPatchScalarField
        (
            const pressureBCFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureBCFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureBCFvPatchScalarField
        (
            const pressureBCFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureBCFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

            const scalarField& freestreamValue() const
            {
                return refValue();
            }

            scalarField& freestreamValue()
            {
                return refValue();
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
