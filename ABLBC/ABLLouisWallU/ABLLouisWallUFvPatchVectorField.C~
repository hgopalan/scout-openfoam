/*---------------------------------------------------------------------------*\
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

#include "ABLLouisWallUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interpolateXY.H"
#include "scalarIOList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  namespace incompressible
  {

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    ABLLouisWallUFvPatchVectorField::
    ABLLouisWallUFvPatchVectorField
    (
     const fvPatch& p,
     const DimensionedField<vector, volMesh>& iF
     )
      :
      fixedValueFvPatchVectorField(p, iF),
      kappa_(0.41),
      z0_(0.1)

    {}


    ABLLouisWallUFvPatchVectorField::
    ABLLouisWallUFvPatchVectorField
    (
     const fvPatch& p,
     const DimensionedField<vector, volMesh>& iF,
     const dictionary& dict
     )
      :
      fixedValueFvPatchVectorField(p, iF),
      kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
      z0_(dict.lookupOrDefault<scalar>("z0", 0.1))
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


    ABLLouisWallUFvPatchVectorField::
    ABLLouisWallUFvPatchVectorField
    (
     const ABLLouisWallUFvPatchVectorField& pvf,
     const fvPatch& p,
     const DimensionedField<vector, volMesh>& iF,
     const fvPatchFieldMapper& mapper
     )
      :
      fixedValueFvPatchVectorField(pvf, p, iF, mapper),
      kappa_(pvf.kappa_),
      z0_(pvf.z0_)
    {}


    ABLLouisWallUFvPatchVectorField::
    ABLLouisWallUFvPatchVectorField
    (
     const ABLLouisWallUFvPatchVectorField& pvf,
     const DimensionedField<vector, volMesh>& iF
     )
      :
      fixedValueFvPatchVectorField(pvf, iF),
      kappa_(pvf.kappa_),
      z0_(pvf.z0_)
    {}


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void ABLLouisWallUFvPatchVectorField::autoMap
    (
     const fvPatchFieldMapper& m
     )
    {
      fixedValueFvPatchVectorField::autoMap(m);

    }


    void ABLLouisWallUFvPatchVectorField::rmap
    (
     const fvPatchVectorField& pvf,
     const labelList& addr
     )
    {
      return;
    }
    void ABLLouisWallUFvPatchVectorField::updateCoeffs()
    {
      if (this->updated())
	{
	  return;
	}
      const scalarIOList& Tvalues=this->db().lookupObject<scalarIOList>("timeProfile");      
      const scalarIOList& MOvalues=this->db().lookupObject<scalarIOList>("MOProfile");      
      const volScalarField p_rgh_=this->db().lookupObject<volScalarField>("p_rgh");  
      scalarField timevalues(Tvalues.size(),0.0);
      scalarField movalues(Tvalues.size(),0.0);      
      forAll(Tvalues,i)
	{
	  timevalues[i] = Tvalues[i];
	  movalues[i] = MOvalues[i];
	}
      scalar MOL=interpolateXY(p_rgh_.time().value(),timevalues,movalues);
    if(mag(MOL)>25)
      {
	MOL=MOL;
      }
    else if(neg(MOL) && MOL> -25)
      {
	MOL=-25;
      }
    else if(pos(MOL) && mag(MOL)<25)
      {
	MOL=25;
      }
      
      // Compute ustar from Louis Formula
      const vectorField  Uc = this->patchInternalField();
      vectorField normalV_=patch().nf();
      scalarField UpParallel=mag(Uc - (Uc & normalV_)*normalV_);
      //const scalarField nhat_=max(mag((patch().Cn()-patch().Cf()) mag((patch().Cn()-patch().Cf()) & patch().nf()) patch().nf()),0.25+0*z_);
      scalarField y_=mag((patch().Cn()-patch().Cf()) & patch().nf());
      const scalarField nhat_=max(y_,0.25+0*y_)+z0_;      
      const volScalarField T_=this->db().lookupObject<volScalarField>("T");
      label patchi=this->patch().index();
      const scalarField Tcell=T_.boundaryField()[patchi].patchInternalField();
      const scalarField TFace=this->patch().lookupPatchField<volScalarField, scalar>("T");
      scalarField deltaTheta=Tcell-TFace;
      scalarField louisA=kappa_/log((nhat_+z0_)/z0_);  //Eq.13
      scalar louisCstar=7.4;
      scalarField Fm=louisF(UpParallel,deltaTheta,Tcell,nhat_,MOL,louisCstar);
      scalarField ustar=louisA*UpParallel*sqrt(Fm); //Eq. 12a
      louisCstar=5.3;
      scalarField Fh=louisF(UpParallel,deltaTheta,Tcell,nhat_,MOL,louisCstar)+SMALL;
      scalarField louisMOL=Tcell/(kappa_*9.81)*ustar*UpParallel/(deltaTheta+SMALL)*Fm/Fh*0.7;
      forAll(louisMOL,faceI)
	{
	  if(mag(louisMOL[faceI]<5.0))
	    {
	      louisMOL[faceI]=MOL;
	    }
	}
      vectorField::operator=(ULouisWall(this->patchInternalField(),nhat_,patch().nf(),ustar,louisMOL));      
      fixedValueFvPatchVectorField::updateCoeffs();
    }

    tmp<vectorField> ABLLouisWallUFvPatchVectorField::ULouisWall(const vectorField& Up,const scalarField& nhat_,const vectorField& nhatV_,scalarField& ustar,scalarField& MOL_) const
    {
      vectorField UpParallel=Up - (Up & nhatV_)*nhatV_;    
      tmp<scalarField> deltaU=ustar/kappa_*nhat_/(nhat_+z0_)*phiM(nhat_,MOL_);
      scalarField blankU(this->size(),1);
      if(this->db().foundObject<volScalarField>("buildingBlanking"))
      {
        const volScalarField building=this->db().lookupObject<volScalarField>("buildingBlanking");
	label patchi=this->patch().index();
        const scalarField buildingBlanking=building.boundaryField()[patchi].patchInternalField();
	blankU=buildingBlanking;
      }
      scalarField blankTree(this->size(),-1);
      if(this->db().foundObject<volScalarField>("treeBlanking"))
	{
	  const volScalarField tree=this->db().lookupObject<volScalarField>("treeBlanking");
	  label patchi=this->patch().index();
	  const scalarField treeBlanking=tree.boundaryField()[patchi].patchInternalField();
	  blankTree=treeBlanking;
	  /*
	  scalar treecounter=0;
	  forAll(blankTree,faceI)
	    {
	      if(blankTree[faceI]>0)
		treecounter=treecounter+1;
	    }
	  reduce(treecounter,sumOp<scalar>());
	  Info<<"Patch:"<<patch().name()<<"  Tree blank:"<<treecounter<<endl;*/
	}
      
      return neg(blankTree)*blankU*UpParallel*(1.0-deltaU/(mag(UpParallel)+SMALL));
    }    
      
    tmp<scalarField> ABLLouisWallUFvPatchVectorField::louisF(scalarField& UpParallel,scalarField& deltaTheta,const scalarField& Tcell,const scalarField& nhat_,scalar& MOL, scalar& louisCstar) const
    {
      // Unstable
      scalarField Rib=9.81*nhat_*deltaTheta/(Tcell*sqr(UpParallel)+SMALL);
      Rib=max(Rib,-10+0*Rib);
      scalar louisB=9.4;
      scalar kappa_=0.41;
      scalarField louisC=louisCstar*louisB*sqr(kappa_/log((nhat_+z0_)/z0_))*sqrt(nhat_/z0_);
      scalar louisBstar=4.7;    
      return neg(Rib)*(1-louisB*Rib/(1+louisC*sqrt(mag(Rib))))+pos(Rib)/sqr(1+louisBstar*Rib+SMALL);
      /*
      // Unstable
      if(MOL<0)
      {
      Rib=neg(Rib)*Rib;
      scalar louisB=9.4;
      scalarField louisC=louisCstar*louisB*sqr(kappa_/log((nhat_+z0_)/z0_))*sqrt(nhat_/z0_);
      return 1-louisB*Rib/(1+louisC*sqrt(mag(Rib)));
      }
      else
      {
      Rib=pos(Rib)*Rib;
      scalar louisBstar=4.7;
      return 1/sqr(1+louisBstar*Rib);
      }*/
    }
    tmp<scalarField> ABLLouisWallUFvPatchVectorField::phiM(const scalarField &p,scalarField& MOL_) const 
    {
      //      return (pos(MOL_)*(1+ 5.0*p/MOL_)+neg(MOL_)*(pow(1- neg(MOL_)*16*p/MOL_,-0.25)));
      return (pos(MOL_)*(1+ 5.0*p/MOL_)+neg(MOL_)*(pow(1- neg(MOL_)*12*pos(p)*p/MOL_,-1.0/3.0)));

    }

    void ABLLouisWallUFvPatchVectorField::write(Ostream& os) const
    {
      fvPatchVectorField::write(os);
      os.writeKeyword("kappa")<< kappa_ << token::END_STATEMENT << nl;
      os.writeKeyword("z0")<< z0_ << token::END_STATEMENT << nl;
      writeEntry("value", os);
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField
    (
     fvPatchVectorField,
     ABLLouisWallUFvPatchVectorField
     );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  }
} // End namespace Foam

// ************************************************************************* //
