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

#include "unsteadyABLInflowURainFvPatchVectorField.H"
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

  unsteadyABLInflowURainFvPatchVectorField::
  unsteadyABLInflowURainFvPatchVectorField
  (
   const fvPatch& p,
   const DimensionedField<vector, volMesh>& iF
   )
    :
    fixedValueFvPatchVectorField(p, iF),
    Vt_(0.0)
  {}

  unsteadyABLInflowURainFvPatchVectorField::
  unsteadyABLInflowURainFvPatchVectorField
  (
   const unsteadyABLInflowURainFvPatchVectorField& pvf,
   const fvPatch& p,
   const DimensionedField<vector, volMesh>& iF,
   const fvPatchFieldMapper& mapper
   )
    :
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    Vt_(pvf.Vt_)
  {}

  unsteadyABLInflowURainFvPatchVectorField::
  unsteadyABLInflowURainFvPatchVectorField
  (
   const fvPatch& p,
   const DimensionedField<vector, volMesh>& iF,
   const dictionary& dict
   )
    :
    fixedValueFvPatchVectorField(p, iF),
    Vt_(dict.lookupOrDefault<scalar>("Vt", 0.0))
  {
    if (dict.found("value"))
      {
	fvPatchField<vector>::operator=
	  (
	   vectorField("value", dict, p.size())
	   );
         
      }
    else
      {
	// Evaluate the wall velocity
	updateCoeffs();
      }
  }




  unsteadyABLInflowURainFvPatchVectorField::
  unsteadyABLInflowURainFvPatchVectorField
  (
   const unsteadyABLInflowURainFvPatchVectorField& pvf
   )
    :
    fixedValueFvPatchVectorField(pvf),
    Vt_(pvf.Vt_)
  {}

  unsteadyABLInflowURainFvPatchVectorField::
  unsteadyABLInflowURainFvPatchVectorField
  (
   const unsteadyABLInflowURainFvPatchVectorField& pvf,
   const DimensionedField<vector, volMesh>& iF
   )
    :
    fixedValueFvPatchVectorField(pvf, iF),
    Vt_(pvf.Vt_)
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  void unsteadyABLInflowURainFvPatchVectorField::autoMap
  (
   const fvPatchFieldMapper& m
   )
  {
    fixedValueFvPatchVectorField::autoMap(m);
  }


  void unsteadyABLInflowURainFvPatchVectorField::rmap
  (
   const fvPatchVectorField& ptf,
   const labelList& addr
   )
  {
    fixedValueFvPatchVectorField::rmap(ptf, addr);
    //      const unsteadyABLInflowURainFvPatchVectorField& nrwfpsf =
    //	refCast<const unsteadyABLInflowURainFvPatchVectorField>(ptf);
  }

  void unsteadyABLInflowURainFvPatchVectorField::updateCoeffs()
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
    fvPatchVectorField::operator=(TRef(flowDir));  
    fixedValueFvPatchVectorField::updateCoeffs();
  }
  tmp<vectorField> unsteadyABLInflowURainFvPatchVectorField::TRef(const vector&flowdir) const
  {
    vectorField patchNormal=patch().nf();
    const scalarField nhat_=patch().Cf().component(vector::Z);
    scalar fluxdir=sum(patchNormal & flowdir);
    if(fluxdir>0)
      {
        return this->patchInternalField();
      }
    const scalarIOList&  zvalues=this->db().lookupObject<scalarIOList>("zProfile");
    const scalarIOList&  Uvalues=this->db().lookupObject<scalarIOList>("vertUProfile");
    const scalarIOList&  Vvalues=this->db().lookupObject<scalarIOList>("vertVProfile");
    scalarField verticalZList(zvalues.size(),0.0);
    scalarField verticalUList(Uvalues.size(),0.0);
    scalarField verticalVList(Vvalues.size(),0.0);
    vector xDir(1,0,0);
    vector yDir(0,1,0);
    vector zDir(0,0,1);
    forAll(zvalues,i)
      {
	verticalZList[i] = zvalues[i];
	verticalUList[i] = Uvalues[i]; 
	verticalVList[i] = Vvalues[i]; 
      }
    scalarField UEqn=0*nhat_+verticalUList[0];
    scalarField VEqn=0*nhat_+verticalVList[0];
    forAll(UEqn,facei)
      {
	UEqn[facei]=interpolateXY(nhat_[facei],verticalZList,verticalUList);
	VEqn[facei]=interpolateXY(nhat_[facei],verticalZList,verticalVList);
	//Info<<patch().name()<<TEqn[facei]<<endl;
      }
    return(xDir*UEqn+yDir*VEqn+zDir*Vt_);

  } 

  void unsteadyABLInflowURainFvPatchVectorField::write(Ostream& os) const
  {
    fvPatchVectorField::write(os);
    os.writeKeyword("Vt")<< Vt_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
  }


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  makePatchTypeField
  (
   fvPatchVectorField,
   unsteadyABLInflowURainFvPatchVectorField
   );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
