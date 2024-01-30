/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "porousMediaAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(porousMediaAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            porousMediaAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::porousMediaAbsorptionEmission::porousMediaAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    a_("absorptivity", coeffsDict_),
    e_("emissivity", coeffsDict_),
    E_("E", coeffsDict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::porousMediaAbsorptionEmission::aCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            a_
        )
    );
    const objectRegistry& db = mesh().thisDb();
    if(db.foundObject<volScalarField>("treeShading") && db.foundObject<volScalarField>("treeA"))
      {
        const volScalarField& treeShading_= db.lookupObject<volScalarField>("treeShading");
        const volScalarField& treeA_= db.lookupObject<volScalarField>("treeA");
        const volScalarField& treeTau_= db.lookupObject<volScalarField>("treeTau");
        return treeA_-treeTau_;
      }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::porousMediaAbsorptionEmission::eCont(const label bandI) const
{
    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            e_
        )
    );
    const objectRegistry& db = mesh().thisDb();
    if(db.foundObject<volScalarField>("treeShading") && db.foundObject<volScalarField>("treeE"))
      {
        const volScalarField& treeShading_= db.lookupObject<volScalarField>("treeShading");
        const volScalarField& treeE_= db.lookupObject<volScalarField>("treeE");
        return treeE_;
      }

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::porousMediaAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            E_
        )
    );
    return tE;
}


// ************************************************************************* //
