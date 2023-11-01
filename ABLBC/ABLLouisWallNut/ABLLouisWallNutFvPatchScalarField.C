/*---------------------------------------------------------------------------*\
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

#include "ABLLouisWallNutFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"
#include "scalarIOList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  ABLLouisWallNutFvPatchScalarField::
  ABLLouisWallNutFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(p, iF),
    Prt_(1.0),
    z0_(0.1)
  {}



  ABLLouisWallNutFvPatchScalarField::
  ABLLouisWallNutFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const dictionary& dict
   )
    :
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt",1.0)),    
    z0_(dict.lookupOrDefault<scalar>("z0",0.1))    
  {

  }

  ABLLouisWallNutFvPatchScalarField::
  ABLLouisWallNutFvPatchScalarField
  (
   const ABLLouisWallNutFvPatchScalarField& ptf,
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const fvPatchFieldMapper& mapper
   )
    :
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
    z0_(ptf.z0_)
  {}


  ABLLouisWallNutFvPatchScalarField::
  ABLLouisWallNutFvPatchScalarField
  (
   const ABLLouisWallNutFvPatchScalarField& rwfpsf,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(rwfpsf, iF),
    Prt_(rwfpsf.Prt_),    
    z0_(rwfpsf.z0_)
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void ABLLouisWallNutFvPatchScalarField::autoMap
  (
   const fvPatchFieldMapper& m
   )
  {
    fixedValueFvPatchScalarField::autoMap(m);
  }


  void ABLLouisWallNutFvPatchScalarField::rmap
  (
   const fvPatchScalarField& ptf,
   const labelList& addr
   )
  {
    return;
  }

  // * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * // 
  void ABLLouisWallNutFvPatchScalarField::updateCoeffs() 
  {
    if (this->updated())
      {
	return;
      }
    const volVectorField& U=this->db().lookupObject<volVectorField>("U");      
    const label patchi = patch().index(); 
    const vectorField  Uc = U.boundaryField()[patchi].patchInternalField();
    const fvPatch& patch=this->patch();
    vectorField normalV_=patch.nf();
    scalarField UpParallel=mag(Uc - (Uc & normalV_)*normalV_);
    scalarField y_=mag((patch.Cn()-patch.Cf()) & patch.nf());
    const scalarField nhat_=max(y_,0.25+0*y_);    
    //const scalarField nhat_=mag((patch.Cn()-patch.Cf()) & patch.nf());
    const scalarIOList& Tvalues=this->db().lookupObject<scalarIOList>("timeProfile");      
    const scalarIOList& MOvalues=this->db().lookupObject<scalarIOList>("MOProfile");      
    const volScalarField rho_=this->db().lookupObject<volScalarField>("rhok");  
    scalarField timevalues(Tvalues.size(),0.0);
    scalarField movalues(Tvalues.size(),0.0);      
    forAll(Tvalues,i)
      {
	timevalues[i] = Tvalues[i];
	movalues[i] = MOvalues[i];
      }
    scalar MOL=interpolateXY(rho_.time().value(),timevalues,movalues);
    if(mag(MOL)>25)
      MOL=MOL;
    else if(neg(MOL) && MOL> -25)
      MOL=-25;
    else if(pos(MOL) && mag(MOL)<25)
      MOL=25;
    // Compute ustar from Louis Formula
    scalar kappa_=0.41;    
    const volScalarField T_=this->db().lookupObject<volScalarField>("T");
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

    //    scalarField ustar=louisA*UpParallel*sqrt(louisF(UpParallel,deltaTheta,Tcell,nhat_,MOL,louisCstar)); //Eq. 12a
    //    louisCstar=5.3;
    //    scalar louisCstarU=7.4;
    //    scalarField thetastar=louisA/0.74*deltaTheta*louisF(UpParallel,deltaTheta,Tcell,nhat_,MOL,louisCstar)/(sqrt(louisF(UpParallel,deltaTheta,Tcell,nhat_,MOL,louisCstarU))+SMALL); //Eq. 12b; //Eq. 12b      
    //    scalarField thetastar=1/max(ustar,0.01+0*ustar)*sqr(louisA)/0.74*UpParallel*deltaTheta*louisF(UpParallel,deltaTheta,Tcell,nhat_,MOL,louisCstar); //Eq. 12b
    //    scalarField louisMOL= MOL+0*Tcell;//sqr(ustar)/(kappa_*9.81*thetastar+SMALL)+SMALL; //Eq.5
    //    scalarField louisMOL= Tcell*sqr(ustar)/(kappa_*9.81*thetastar+SMALL)+SMALL; //Eq.5
    //    scalarField nutw=ustar*kappa_*(nhat_+z0_)/(phiM(nhat_,louisMOL)+SMALL)*1.0/Prt_;
    //        const scalarField blockOff=this->patch().lookupPatchField<volScalarField, scalar>("unphysicalCells");
    scalarField nutw=ustar*kappa_*z0_/(phiM(z0_+0*nhat_,louisMOL)+SMALL)*1.0/Prt_;
    scalarField blankNut(this->size(),1.0);
    if(this->db().foundObject<volScalarField>("buildingBlanking"))
      {
        const volScalarField building=this->db().lookupObject<volScalarField>("buildingBlanking");
	label patchi=this->patch().index();
        const scalarField buildingBlanking=building.boundaryField()[patchi].patchInternalField();
        blankNut=buildingBlanking;
      }
    scalarField blankTree(this->size(),-1);
    if(this->db().foundObject<volScalarField>("treeBlanking"))
      {
	const volScalarField tree=this->db().lookupObject<volScalarField>("treeBlanking");
	label patchi=this->patch().index();
	const scalarField treeBlanking=tree.boundaryField()[patchi].patchInternalField();
	blankTree=treeBlanking;
      }
    scalarField::operator=(neg(blankTree)*(blankNut*nutw+(1-blankNut)*1e-5));
    fixedValueFvPatchScalarField::updateCoeffs();
  }

  tmp<scalarField> ABLLouisWallNutFvPatchScalarField::louisF(scalarField& UpParallel,scalarField& deltaTheta,const scalarField& Tcell,const scalarField& nhat_,scalar& MOL, scalar& louisCstar) const
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
    scalar kappa_=0.41;
    scalarField louisC=louisCstar*louisB*sqr(kappa_/log((nhat_+z0_)/z0_))*sqrt(nhat_/z0_);
    return 1-louisB*Rib/(1+louisC*sqrt(mag(Rib)));
    }
    else
    {
    Rib=pos(Rib)*Rib;
    scalar louisBstar=4.7;
    return 1/sqr(1+louisBstar*Rib);
    } */
  }
  // Stability Function 
  tmp<scalarField> ABLLouisWallNutFvPatchScalarField::phiM(const scalarField &p,scalarField& MOL_) const 
  {
    //      return (pos(MOL_)*(1+ 5.0*p/MOL_)+neg(MOL_)*(pow(1- neg(MOL_)*16*pos(p)*p/MOL_,-0.25)));
    return (pos(MOL_)*(1+ 5.0*p/MOL_)+neg(MOL_)*(pow(1- neg(MOL_)*12*pos(p)*p/MOL_,-1.0/3.0)));
  }

  void ABLLouisWallNutFvPatchScalarField::write(Ostream& os) const
  {
    //      fvPatchField<scalar>::write(os);
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("Prt")<< Prt_ << token::END_STATEMENT << nl;    
    os.writeKeyword("z0")<< z0_ << token::END_STATEMENT << nl;    
    //      writeEntry("value", os);
  }

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  makePatchTypeField
  (
   fvPatchScalarField,
   ABLLouisWallNutFvPatchScalarField
   );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
