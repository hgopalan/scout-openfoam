/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "ABLMOWallQvFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"
#include "scalarIOList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  ABLMOWallQvFvPatchScalarField::
  ABLMOWallQvFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(p, iF),
    z0_(0.1),
    mavail(0.1)
  {}



  ABLMOWallQvFvPatchScalarField::
  ABLMOWallQvFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const dictionary& dict
   )
    :
    fixedValueFvPatchScalarField(p, iF, dict),
    z0_(dict.lookupOrDefault<scalar>("z0",0.1)),
    mavail(dict.lookupOrDefault<scalar>("mavail",0.1))          
  {
    fvPatchField<scalar>::operator=(patchInternalField());
    /*
      else
      {
      // Evaluate the wall velocity
      updateCoeffs();
      } */
  }

  ABLMOWallQvFvPatchScalarField::
  ABLMOWallQvFvPatchScalarField
  (
   const ABLMOWallQvFvPatchScalarField& ptf,
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const fvPatchFieldMapper& mapper
   )
    :
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    z0_(ptf.z0_),
    mavail(ptf.mavail)
      
  {}


  ABLMOWallQvFvPatchScalarField::
  ABLMOWallQvFvPatchScalarField
  (
   const ABLMOWallQvFvPatchScalarField& rwfpsf,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(rwfpsf, iF),
    z0_(rwfpsf.z0_),
    mavail(rwfpsf.mavail)
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void ABLMOWallQvFvPatchScalarField::autoMap
  (
   const fvPatchFieldMapper& m
   )
  {
    fixedValueFvPatchScalarField::autoMap(m);
  }


  void ABLMOWallQvFvPatchScalarField::rmap
  (
   const fvPatchScalarField& ptf,
   const labelList& addr
   )
  {
    fixedValueFvPatchScalarField::rmap(ptf, addr);



  }

  // * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * // 
  void ABLMOWallQvFvPatchScalarField::updateCoeffs() 
  {
    if (this->updated())
      {
	return;
      }
    const volScalarField p_rgh_=this->db().lookupObject<volScalarField>("p_rgh"); 
    // The reference wind speed is read from constant/ABLDict. 
    // Default: (Time U V T MOL)
    const scalarIOList&  Tvalues=this->db().lookupObject<scalarIOList>("timeProfile");
    const scalarIOList& Uvalues=this->db().lookupObject<scalarIOList>("uProfile");
    const scalarIOList& Vvalues=this->db().lookupObject<scalarIOList>("vProfile");
    const scalarIOList& MOvalues=this->db().lookupObject<scalarIOList>("MOProfile");
    const scalarIOList& PBLHvalues=this->db().lookupObject<scalarIOList>("PBLHProfile");
    const scalarIOList& MBLvalues=this->db().lookupObject<scalarIOList>("MBLProfile");
    scalarField timevalues(Tvalues.size(),0.0);
    scalarField uvalues(Tvalues.size(),0.0);
    scalarField vvalues(Tvalues.size(),0.0);
    scalarField movalues(Tvalues.size(),0.0);
    scalarField pblhvalues(Tvalues.size(),0.0);
    scalarField mblvalues(Tvalues.size(),0.0);
    forAll(Tvalues,i)
      {
	timevalues[i] = Tvalues[i];
	uvalues[i] = Uvalues[i];
	vvalues[i] = Vvalues[i];
	movalues[i] = MOvalues[i];
	pblhvalues[i] = PBLHvalues[i];
	mblvalues[i] = MBLvalues[i];
      }
    scalar u1=interpolateXY(p_rgh_.time().value(),timevalues,uvalues);
    scalar v1=interpolateXY(p_rgh_.time().value(),timevalues,vvalues);
    scalar M=sqrt(sqr(u1)+sqr(v1));
    scalar MOL=interpolateXY(p_rgh_.time().value(),timevalues,movalues);
    if(mag(MOL)>25)
      MOL=MOL;
    else if(neg(MOL) && MOL> -25)
      MOL=-25;
    else if(pos(MOL) && mag(MOL)<25)
      MOL=25;
    
    scalar PBLH=interpolateXY(p_rgh_.time().value(),timevalues,pblhvalues);
    scalar MBL=interpolateXY(p_rgh_.time().value(),timevalues,mblvalues);
    const vector flowDir(u1/M,v1/M,0);
    scalarField qVw(this->size(),1e-4);
    // Neutral
    if(mag(MOL)>500)
      {
	qVw=(neutral(M,flowDir,MOL,PBLH,MBL)); 
      }
    else if(pos(MOL) && MOL<2.5)
      {
	qVw=(veryStable(M,flowDir,MOL,PBLH,MBL)); 
      }
    else if(pos(MOL))
      {
	qVw=(stable(M,flowDir,MOL,PBLH,MBL)); 
      }
    else 
      {
	qVw=(unstable(M,flowDir,MOL,PBLH,MBL)); 
      }
    scalarField blankqV(this->size(),1.0);
    if(this->db().foundObject<volScalarField>("buildingBlanking"))
      {
        const volScalarField building=this->db().lookupObject<volScalarField>("buildingBlanking");
	label patchi=patch().index();
        const scalarField buildingBlanking=building.boundaryField()[patchi].patchInternalField();
        blankqV=buildingBlanking;
      }   
    scalarField::operator=(blankqV*qVw+(1-blankqV)*1e-4);
    //scalar qMin_=min(qVw);                                                                                  
    //scalar qMax_=max(qVw);
    //reduce(qMin_, minOp<scalar>());
    //reduce(qMax_, maxOp<scalar>());
    //Info<<this->patch().name()<<" "<<":[ "<<qMin_<<"  "<<qMax_<<" ]"<<endl;    
    fixedValueFvPatchScalarField::updateCoeffs();
  }

  tmp<scalarField> ABLMOWallQvFvPatchScalarField::neutral(scalar M,vector flowDir,scalar MOL, scalar PBLH, scalar MBL) const
  {
    scalar kappa_=0.41;
    // Read wind speed at cell center from registry 
    const volVectorField& U=this->db().lookupObject<volVectorField>("U");
    const label patchi = patch().index();
    const vectorField  Uc = U.boundaryField()[patchi].patchInternalField();
    // Compute u* from Cell-Center
    scalarField y_=mag((patch().Cn()-patch().Cf()) & patch().nf());
    scalarField z_=max(y_,0.25+0*y_)+z0_;
    scalarField denominator=log(z_/z0_)+z_/MBL-z_/PBLH*0.5*z_/MBL;
    vectorField UpParallel=Uc - (Uc & patch().nf())*patch().nf();
    scalarField ustar=mag(UpParallel)*kappa_/denominator;
    scalarField Tw_=patch().lookupPatchField<volScalarField, scalar>("T"); 
    scalarField qvCell=this->patchInternalField();
    scalarField faceesat=610.94*exp(17.625*(Tw_-273.15)/(Tw_-30.11));
    scalarField qvFace=(0.622 *faceesat/(101325-0.378*faceesat));
    scalar ka=2.4e-5;
    scalarField den=log(kappa_*ustar*z_/ka+z_/0.01);
    //scalarField qstar=mavail*kappa_/den*(qvCell-qvFace);
    const volScalarField treeTemp=this->db().lookupObject<volScalarField>("treeBlanking");
    const scalarField treeCell=treeTemp.boundaryField()[patchi].patchInternalField();    
    scalarField qstar=(pos(treeCell)+neg(treeCell)*mavail)*kappa_/den*(qvCell-qvFace);
    scalarField faceqV=qvCell-qstar/kappa_;
    if(mavail>0.9)
      faceqV=qvFace;    
    return max(faceqV,qvCell);
  }
  
  tmp<scalarField> ABLMOWallQvFvPatchScalarField::veryStable(scalar M,vector flowDir,scalar MOL, scalar PBLH, scalar MBL) const
  {
    scalar kappa_=0.41;
    // Read wind speed at cell center from registry 
    const volVectorField& U=this->db().lookupObject<volVectorField>("U");
    const label patchi = patch().index();
    const vectorField  Uc = U.boundaryField()[patchi].patchInternalField();
    // Compute u* from Cell-Center
    //scalarField z_=max(mag((patch().Cn()-patch().Cf()) mag((patch().Cn()-patch().Cf()) & patch().nf()) patch().nf()),0.25+0*z_)+z0_;
    scalarField y_=mag((patch().Cn()-patch().Cf()) & patch().nf());
    scalarField z_=max(y_,0.25+0*y_)+z0_;    
    scalarField denominator=6*log(z_/z0_)-5*z_/PBLH+z_/MBL-z_/PBLH*0.5*z_/MBL;
    vectorField UpParallel=Uc - (Uc & patch().nf())*patch().nf();
    scalarField ustar=mag(UpParallel)*kappa_/denominator;
    scalarField Tw_=patch().lookupPatchField<volScalarField, scalar>("T"); 
    scalarField qvCell=this->patchInternalField();
    scalarField faceesat=610.94*exp(17.625*(Tw_-273.15)/(Tw_-30.11));
    scalarField qvFace=(0.622 *faceesat/(101325-0.378*faceesat));
    scalar ka=2.4e-5;
    scalarField den=log(kappa_*ustar*z_/ka+z_/0.01)+5.0;
    //  scalarField qstar=mavail*kappa_/den*(qvCell-qvFace);
    const volScalarField treeTemp=this->db().lookupObject<volScalarField>("treeBlanking");
    const scalarField treeCell=treeTemp.boundaryField()[patchi].patchInternalField();    
    scalarField qstar=(pos(treeCell)+neg(treeCell)*mavail)*kappa_/den*(qvCell-qvFace);    
    scalarField faceqV=qvCell-qstar/kappa_*6;
    if(mavail>0.9)
      faceqV=qvFace;    
    return max(faceqV,qvCell);
  }
  tmp<scalarField> ABLMOWallQvFvPatchScalarField::stable(scalar M,vector flowDir,scalar MOL, scalar PBLH, scalar MBL) const
  {
    scalar kappa_=0.41;
    // Read wind speed at cell center from registry 
    const volVectorField& U=this->db().lookupObject<volVectorField>("U");
    const label patchi = patch().index();
    const vectorField  Uc = U.boundaryField()[patchi].patchInternalField();
    // Compute u* from Cell-Center
    //scalarField z_=max(mag((patch().Cn()-patch().Cf()) mag((patch().Cn()-patch().Cf()) & patch().nf()) patch().nf()),0.25+0*z_)+z0_;
    scalarField y_=mag((patch().Cn()-patch().Cf()) & patch().nf());
    scalarField z_=max(y_,0.25+0*y_)+z0_;    
    scalarField denominator=log(z_/z0_)+5*z_/MOL*(1-0.5*z_/PBLH)+z_/MBL-z_/PBLH*0.5*z_/MBL;
    vectorField UpParallel=Uc - (Uc & patch().nf())*patch().nf();
    scalarField ustar=mag(UpParallel)*kappa_/denominator;
    // Read wind speed at cell center from registry 
    scalarField Tw_=patch().lookupPatchField<volScalarField, scalar>("T"); 
    scalarField qvCell=this->patchInternalField();
    scalarField faceesat=610.94*exp(17.625*(Tw_-273.15)/(Tw_-30.11));
    scalarField qvFace=(0.622 *faceesat/(101325-0.378*faceesat));
    scalar ka=2.4e-5;
    scalarField den=log(kappa_*ustar*z_/ka+z_/0.01)+5*z_/MOL;
    //  scalarField qstar=mavail*kappa_/den*(qvCell-qvFace);
    const volScalarField treeTemp=this->db().lookupObject<volScalarField>("treeBlanking");
    const scalarField treeCell=treeTemp.boundaryField()[patchi].patchInternalField();    
    scalarField qstar=(pos(treeCell)+neg(treeCell)*mavail)*kappa_/den*(qvCell-qvFace);    
    scalarField phiM=1+5*z_/MOL;
    scalarField faceqV=qvCell-qstar/kappa_*phiM;
    if(mavail>0.9)
      faceqV=qvFace;    
    return max(faceqV,qvCell);
  }
  
  tmp<scalarField> ABLMOWallQvFvPatchScalarField::unstable(scalar M,vector flowDir,scalar MOL, scalar PBLH, scalar MBL) const
  {
    scalar kappa_=0.41;
    // Read wind speed at cell center from registry 
    const volVectorField& U=this->db().lookupObject<volVectorField>("U");
    const label patchi = patch().index();
    const vectorField  Uc = U.boundaryField()[patchi].patchInternalField();
    // Compute u* from Cell-Center
    //scalarField z_=max(mag((patch().Cn()-patch().Cf()) & patch().nf()),0.25+0*z_)+z0_;
    scalarField y_=mag((patch().Cn()-patch().Cf()) & patch().nf());
    scalarField z_=max(y_,0.25+0*y_)+z0_;    
    //    Pout<<min(z_)<<"  "<<max(z_)<<endl;
    scalarField xF=pow(1-12*z_/MOL,1.0/3.0);
    scalarField psiMF=1.5*log((1+xF+sqr(xF))/3.0)+3.14/sqrt(3.0)-sqrt(3.0)*atan((1+2*xF)/sqrt(3.0));
    //    Pout<<"psiMF:"<<psiMF<<endl;
    scalarField denominator=log(z_/z0_)-psiMF+z_/PBLH*(1+(pow(1-12*z_/MOL,1.0/3.0)-1)/(8*z_/MOL));
    denominator=denominator+z_/MBL-z_/PBLH*0.5*z_/MBL;
    vectorField UpParallel=Uc- (Uc & patch().nf())*patch().nf();
    scalarField ustar=mag(UpParallel)*kappa_/denominator;
    //    Pout<<"denominator:"<<psiMF<<endl;
    // Read wind speed at cell center from registry 
    scalarField Tw_=patch().lookupPatchField<volScalarField, scalar>("T");   
    scalarField qvCell=this->patchInternalField();
    scalarField faceesat=610.94*exp(17.625*(Tw_-273.15)/(Tw_-30.11));
    scalarField qvFace=(0.622 *faceesat/(101325-0.378*faceesat));
    scalar ka=2.4e-5;
    scalarField x=pow(1- 16*z_/MOL,0.25);
    scalarField den=log(kappa_*ustar*z_/ka+z_/0.01)-2*log(0.5*(1+sqr(x)));
    //  scalarField qstar=mavail*kappa_/den*(qvCell-qvFace);
    const volScalarField treeTemp=this->db().lookupObject<volScalarField>("treeBlanking");
    const scalarField treeCell=treeTemp.boundaryField()[patchi].patchInternalField();    
    scalarField qstar=(pos(treeCell)+neg(treeCell)*mavail)*kappa_/den*(qvCell-qvFace);    
    scalarField phiM=pow(1- 16*z_/MOL,-1.0/2.0);
    scalarField faceqV=qvCell-qstar/kappa_*phiM;
    if(mavail>0.9)
      faceqV=qvFace;
    return max(faceqV,qvCell);
  }  
  void ABLMOWallQvFvPatchScalarField::write(Ostream& os) const
  {
    fvPatchField<scalar>::write(os);
    os.writeKeyword("z0")<< z0_ << token::END_STATEMENT << nl;
    os.writeKeyword("mavail")<< mavail << token::END_STATEMENT << nl;    
    writeEntry("value", os);
  }


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  makePatchTypeField
  (
   fvPatchScalarField,
   ABLMOWallQvFvPatchScalarField
   );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// ************************************************************************* //
