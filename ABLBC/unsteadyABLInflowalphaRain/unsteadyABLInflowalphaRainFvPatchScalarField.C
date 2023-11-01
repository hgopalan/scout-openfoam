/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
  \\/     M anipulation  |
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

#include "unsteadyABLInflowalphaRainFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interpolateXY.H"
#include "scalarIOList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  unsteadyABLInflowalphaRainFvPatchScalarField::
  unsteadyABLInflowalphaRainFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(p, iF),
    alphaValue(0.0)
  {}

  unsteadyABLInflowalphaRainFvPatchScalarField::
  unsteadyABLInflowalphaRainFvPatchScalarField
  (
   const unsteadyABLInflowalphaRainFvPatchScalarField& pvf,
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const fvPatchFieldMapper& mapper
   )
    :
    fixedValueFvPatchScalarField(pvf, p, iF, mapper),
    alphaValue(pvf.alphaValue)
  {}

  unsteadyABLInflowalphaRainFvPatchScalarField::
  unsteadyABLInflowalphaRainFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const dictionary& dict
   )
    :
    fixedValueFvPatchScalarField(p, iF),
    alphaValue(dict.lookupOrDefault<scalar>("alphaValue", 0.0))
  {
    if (dict.found("value"))
      {
	fvPatchField<scalar>::operator=
	  (
	   scalarField("value", dict, p.size())
	   );
         
      }
    else
      {
	// Evaluate the wall velocity
	updateCoeffs();
      }
  }




  unsteadyABLInflowalphaRainFvPatchScalarField::
  unsteadyABLInflowalphaRainFvPatchScalarField
  (
   const unsteadyABLInflowalphaRainFvPatchScalarField& pvf
   )
    :
    fixedValueFvPatchScalarField(pvf),
    alphaValue(pvf.alphaValue)
  {}

  unsteadyABLInflowalphaRainFvPatchScalarField::
  unsteadyABLInflowalphaRainFvPatchScalarField
  (
   const unsteadyABLInflowalphaRainFvPatchScalarField& pvf,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(pvf, iF),
    alphaValue(pvf.alphaValue)
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  void unsteadyABLInflowalphaRainFvPatchScalarField::autoMap
  (
   const fvPatchFieldMapper& m
   )
  {
    fixedValueFvPatchScalarField::autoMap(m);
  }


  void unsteadyABLInflowalphaRainFvPatchScalarField::rmap
  (
   const fvPatchScalarField& ptf,
   const labelList& addr
   )
  {
    fixedValueFvPatchScalarField::rmap(ptf, addr);
    //      const unsteadyABLInflowalphaRainFvPatchScalarField& nrwfpsf =
    //	refCast<const unsteadyABLInflowalphaRainFvPatchScalarField>(ptf);
  }

  void unsteadyABLInflowalphaRainFvPatchScalarField::updateCoeffs()
  {
    if (this->updated())
      {
	return;
      }
    const volScalarField p_rgh_=this->db().lookupObject<volScalarField>("p_rgh");  
    const scalarIOList&  Tvalues=this->db().lookupObject<scalarIOList>("timeProfile");
    const scalarIOList& Uvalues=this->db().lookupObject<scalarIOList>("uProfile");
    const scalarIOList& Vvalues=this->db().lookupObject<scalarIOList>("vProfile");
    scalarField timevalues(Tvalues.size(),0.0);
    scalarField uvalues(Tvalues.size(),0.0);
    scalarField vvalues(Tvalues.size(),0.0);
    forAll(Tvalues,i)
      {
	timevalues[i] = Tvalues[i];
	uvalues[i] = Uvalues[i];
	vvalues[i] = Vvalues[i];
      }
    scalar u1=interpolateXY(p_rgh_.time().value(),timevalues,uvalues);
    scalar v1=interpolateXY(p_rgh_.time().value(),timevalues,vvalues);
    scalar M=sqrt(sqr(u1)+sqr(v1));
    const vector flowDir(u1/M,v1/M,0);
    fvPatchScalarField::operator=(TRef(flowDir));  
    fixedValueFvPatchScalarField::updateCoeffs();
  }
  tmp<scalarField> unsteadyABLInflowalphaRainFvPatchScalarField::TRef(const vector&flowdir) const
  {
    vectorField patchNormal=patch().nf();
    const scalarField nhat_=patch().Cf().component(vector::Z);
    scalar fluxdir=sum(patchNormal & flowdir);
    if(fluxdir>0)
      {
        return this->patchInternalField();
      }
    return(0*this->patchInternalField()+alphaValue);

  } 

  void unsteadyABLInflowalphaRainFvPatchScalarField::write(Ostream& os) const
  {
    fvPatchScalarField::write(os);
    os.writeKeyword("alphaValue")<< alphaValue << token::END_STATEMENT << nl;
    writeEntry("value", os);
  }


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  makePatchTypeField
  (
   fvPatchScalarField,
   unsteadyABLInflowalphaRainFvPatchScalarField
   );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
