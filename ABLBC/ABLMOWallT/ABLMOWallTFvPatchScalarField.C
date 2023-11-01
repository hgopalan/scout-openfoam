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

#include "ABLMOWallTFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"
#include "scalarIOList.H"
#include "pimpleControl.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


 

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  ABLMOWallTFvPatchScalarField::
  ABLMOWallTFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(p, iF),
    z0_(0.1),
    matDensity_(),
    specificHeat_(),
    thermalCon_(),
    thickness_()
  {}


  ABLMOWallTFvPatchScalarField::
  ABLMOWallTFvPatchScalarField
  (
   const ABLMOWallTFvPatchScalarField& ptf,
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const fvPatchFieldMapper& mapper
   )
    :
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    z0_(ptf.z0_),
    matDensity_(ptf.matDensity_),
    specificHeat_(ptf.specificHeat_),
    thermalCon_(ptf.thermalCon_),
    thickness_(ptf.thickness_),
    T1_(),
    T2_(),
    T3_(),
    T4_(),
    Tdeep_()
  {
    T1_.setSize(mapper.size());
    T1_.map(ptf.T1_, mapper);    
    T2_.setSize(mapper.size());
    T2_.map(ptf.T2_, mapper);
    T3_.setSize(mapper.size());
    T3_.map(ptf.T3_, mapper);
    T4_.setSize(mapper.size());
    T4_.map(ptf.T4_, mapper);
    Tdeep_.setSize(mapper.size());
    Tdeep_.map(ptf.Tdeep_, mapper);
  }
  ABLMOWallTFvPatchScalarField::
  ABLMOWallTFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const dictionary& dict
   )
    :
    fixedValueFvPatchScalarField(p, iF),
    z0_(dict.lookupOrDefault<scalar>("z0",0.1)),
    matDensity_(),
    specificHeat_(),
    thermalCon_(),
    thickness_()
   
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
        fvPatchField<scalar>::operator=(patchInternalField());
	//	fvPatchField<scalar>::operator=(this->patchInternalField());
	//	updateCoeffs();
      }
    dict.readEntry("thermalCon",thermalCon_);
    dict.readEntry("thickness",thickness_);
    dict.readEntry("matDensity",matDensity_);
    dict.readEntry("specificHeat",specificHeat_);
    if (dict.found("T1"))
      {
        T1_= scalarField("T1", dict, p.size());
      }
    if (dict.found("T2"))
      {
        T2_= scalarField("T2", dict, p.size());
      }
    if (dict.found("T3"))
      {
        T3_= scalarField("T3", dict, p.size());
      }
    if (dict.found("T4"))
      {
        T4_= scalarField("T4", dict, p.size());
      }
    if (dict.found("Tdeep"))
      {
        Tdeep_= scalarField("Tdeep", dict, p.size());
      }
    

  }


  ABLMOWallTFvPatchScalarField::
  ABLMOWallTFvPatchScalarField
  (
   const ABLMOWallTFvPatchScalarField& rwfpsf
   )
    :
    fixedValueFvPatchScalarField(rwfpsf),
    z0_(rwfpsf.z0_),
    matDensity_(rwfpsf.matDensity_),
    specificHeat_(rwfpsf.specificHeat_),
    thermalCon_(rwfpsf.thermalCon_),
    thickness_(rwfpsf.thickness_),
    T1_(rwfpsf.T1_),
    T2_(rwfpsf.T2_),
    T3_(rwfpsf.T3_),
    T4_(rwfpsf.T4_),
    Tdeep_(rwfpsf.Tdeep_)   
  {}



  ABLMOWallTFvPatchScalarField::
  ABLMOWallTFvPatchScalarField
  (
   const ABLMOWallTFvPatchScalarField& rwfpsf,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(rwfpsf, iF),
    z0_(rwfpsf.z0_),
    matDensity_(rwfpsf.matDensity_),
    specificHeat_(rwfpsf.specificHeat_),    
    thermalCon_(rwfpsf.thermalCon_),
    thickness_(rwfpsf.thickness_),
    T1_(rwfpsf.T1_),
    T2_(rwfpsf.T2_),
    T3_(rwfpsf.T3_),
    T4_(rwfpsf.T4_),
    Tdeep_(rwfpsf.Tdeep_)    

  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void ABLMOWallTFvPatchScalarField::autoMap
  (
   const fvPatchFieldMapper& m
   )
  {
    fixedValueFvPatchScalarField::autoMap(m);
    T1_.autoMap(m);
    T2_.autoMap(m);
    T3_.autoMap(m);
    T4_.autoMap(m);
    Tdeep_.autoMap(m);
  }


  void ABLMOWallTFvPatchScalarField::rmap
  (
   const fvPatchScalarField& ptf,
   const labelList& addr
   )
  {
    fixedValueFvPatchScalarField::rmap(ptf, addr);
    const ABLMOWallTFvPatchScalarField& ewhftpsf =
      refCast<const ABLMOWallTFvPatchScalarField>(ptf);
    T1_.rmap(ewhftpsf.T1_, addr);
    T2_.rmap(ewhftpsf.T2_, addr);
    T3_.rmap(ewhftpsf.T3_, addr);
    T4_.rmap(ewhftpsf.T4_, addr);
    Tdeep_.rmap(ewhftpsf.Tdeep_, addr);
    return;
  }

  // * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

  void ABLMOWallTFvPatchScalarField::updateCoeffs() 
  { 
    if (this->updated())
      {
	return;
      }
    const fvPatch& patch=this->patch();
    /*// Short-Wave
      scalarField qSW(this->size(), 0.0);
      // Short-Wave is stored in QrSolar. 
      if(this->db().foundObject<volScalarField>("qrSolar"))
      {      
      qSW = patch.lookupPatchField<volScalarField, scalar>("qrSolar");
      }    
      // Long-Wave (Built-Environment)
      scalarField qLW(this->size(), 0.0);
      if(this->db().foundObject<volScalarField>("qr"))
      {      
      qLW = patch.lookupPatchField<volScalarField, scalar>("qr");
      }
      if(this->db().foundObject<volScalarField>("Qr"))
      {      
      qLW = patch.lookupPatchField<volScalarField, scalar>("Qr");
      }
      // Long-Wave (sky)
      vector verticalDir(0,0,1);
      scalarField cosEpsilon=(verticalDir & -patch.nf());
      // Incoming Radiation from Clouds (not accounted in the exchanges)
      //    scalarField Tw_=patch.lookupPatchField<volScalarField, scalar>("T");
      scalarField qSky_=0.5*(1+cosEpsilon)*(213+5.5*(this->patchInternalField()-273.15));*/
    scalarField radiation_=(totalRadiation());
    // Sensible Heat Flux 
    const label patchi = patch.index();     
    const volVectorField& U=this->db().lookupObject<volVectorField>("U");
    const vectorField  Uc = U.boundaryField()[patchi].patchInternalField();    
    vectorField normalV_=patch.nf();
    //scalarField nhat_=mag((patch.Cn()-patch.Cf()) & patch.nf());
    scalarField y_=mag((patch.Cn()-patch.Cf()) & patch.nf());
    scalarField nhat_=max(y_,0.25+0*y_);//+z0_;    
    const scalarIOList& Tvalues=this->db().lookupObject<scalarIOList>("timeProfile");      
    const scalarIOList& MOvalues=this->db().lookupObject<scalarIOList>("MOProfile");
    const scalarIOList& PBLHvalues=this->db().lookupObject<scalarIOList>("PBLHProfile");          
    const volScalarField p_rgh_=this->db().lookupObject<volScalarField>("p_rgh");  
    scalarField timevalues(Tvalues.size(),0.0);
    scalarField movalues(Tvalues.size(),0.0);
    scalarField pblhvalues(Tvalues.size(),0.0);
    scalarField switchprofilevalues(Tvalues.size(),0.0);
    forAll(Tvalues,i)
      {
        timevalues[i] = Tvalues[i];
        movalues[i] = MOvalues[i];
        pblhvalues[i] = PBLHvalues[i];	
      }
    scalar kappa_=0.41;
    scalar MOL=interpolateXY(p_rgh_.time().value(),timevalues,movalues);
    if(mag(MOL)>25)
      MOL=MOL;
    else if(neg(MOL) && MOL> -25)
      MOL=-25;
    else if(pos(MOL) && mag(MOL)<25)
      MOL=25;
    
    scalar PBLH=interpolateXY(p_rgh_.time().value(),timevalues,pblhvalues);
    scalarField psiML_=psiM(nhat_,MOL)+SMALL;
    scalarField phiML_=phiM(nhat_,MOL)+SMALL;
    vectorField UParallel= Uc - (Uc & normalV_)*normalV_; 
    scalarField ustar=mag(UParallel)*kappa_/(log((nhat_+z0_)/z0_)-psiML_+SMALL);      
    scalar C=0.01;
    scalarField wstar=pow(C*9.81/this->patchInternalField()*PBLH*pos(radiation_)*(radiation_)/1225,0.33);
    wstar=sqrt(sqr(ustar)+sqr(wstar));
    wstar=min(wstar,0.1+0*wstar);// Suggestion from WRF 
    wstar=max(wstar,1.0+0*wstar); 
    // Latent Heat Flux
    scalarField latentHeat_(this->size(),0.0);
    if(this->db().foundObject<volScalarField>("qv"))
      {
	latentHeat_=(latentHeat(ustar,nhat_,MOL));
      }
    scalar dt_=1e30;
    if((word(U.mesh().ddtScheme("default")) != "steadyState"))
      dt_=U.mesh().time().deltaTValue();
    scalarField Twall=(this->patchInternalField());
    // First-layer 
    scalarField groundHeat_=(groundHeat());
    // Surface Energy Balance
    // R=SH+LH+G 
    // Incoming is positive and outgoing is negative 
    scalarField q_=radiation_-latentHeat_-groundHeat_;
    // Sensible heat flux is rho cp ustar thetastar
    scalarField thetastar= q_/(1225*wstar);
    // thetastar * kappa / zh * phiH = Tw - Tc /zh
    Twall=this->patchInternalField()+thetastar/kappa_*nhat_/(nhat_+z0_)*phiH(nhat_,MOL);
    // Wall heat conduction
    T1_=firstLayer(Twall,dt_);
    // Second-layer
    T2_=secondLayer(dt_);
    // Third-layer
    T3_=thirdLayer(dt_);
    // Last-Layer
    T4_=fourthLayer(dt_);
    /*scalar qMin_=min(Twall);
      scalar qMax_=max(Twall);
      reduce(qMin_, minOp<scalar>());
      reduce(qMax_, maxOp<scalar>());
      Info<<patch.name()<<" Tw  "<<qMin_<<"  "<<qMax_<<endl;
      qMin_=min(ustar);
      qMax_=max(ustar);
      reduce(qMin_, minOp<scalar>());
      reduce(qMax_, maxOp<scalar>());
      Info<<patch.name()<<" w*  "<<qMin_<<"  "<<qMax_<<endl;
      qMin_=min(qLH+qGround_);
      qMax_=max(qLH+qGround_);
      reduce(qMin_, minOp<scalar>());
      reduce(qMax_, maxOp<scalar>());
      Info<<patch.name()<<"Radiation:  "<<qMin_<<"  "<<qMax_<<endl; */
    scalarField::operator=(Twall);
    fixedValueFvPatchScalarField::updateCoeffs();

  }
  // Stability Function 
  tmp<scalarField> ABLMOWallTFvPatchScalarField::psiM(scalarField &p,scalar& MOL_) const
  {
    scalarField x=pow(1-neg(MOL_)*12*p/MOL_,1.0/3.0);
    return (pos(MOL_)*(-5.0* p/MOL_)+neg(MOL_)*(1.5*log(1.0/3.0*(1+x+sqr(x)))-sqrt(3.0)*atan((1+2*x)/sqrt(3.0))+3.14/sqrt(3.0)));

    /*if(pos(MOL_))
      return (-5.0* p/MOL_);
      else
      {
      scalarField x=pow(1- 16*p/MOL_,0.25);
      return (log(0.5*(1+sqr(x))*sqr(0.5*(1+x)))-2*atan(x)+0.5*3.14);
      }*/
  }
    
  tmp<scalarField> ABLMOWallTFvPatchScalarField::psiH(scalarField &p,scalar& MOL_) const
  {
    if (MOL_ < 0)
      {
	scalarField x=pow(1- 16*(p)/MOL_,0.25);
	return 2*log(0.5*(1+sqr(x)));
      }
    else
      return (-5.0*p/MOL_);
  }

  // Stability Function 
  tmp<scalarField> ABLMOWallTFvPatchScalarField::phiM(scalarField &p,scalar& MOL_) const 
  {
    /*if(pos(MOL_))
      return (1+ 5.0*p/MOL_);
      else
      return (pow(1- 16*p/MOL_,-0.25));*/
    return (pos(MOL_)*(1+ 5.0*p/MOL_)+neg(MOL_)*(pow(1- neg(MOL_)*12*pos(p)*p/MOL_,-1.0/3.0)));
  }

  tmp<scalarField> ABLMOWallTFvPatchScalarField::phiH(scalarField &p,scalar& MOL_) const
  {
    if(pos(MOL_))
      return (1+ 5.0*p/MOL_);
    else
      return (pow(1- 16*p/MOL_,-0.5));
  }
  tmp<scalarField> ABLMOWallTFvPatchScalarField::totalRadiation() const
  {
    const fvPatch& patch=this->patch();
    // Short-Wave
    scalarField qSW(this->size(), 0.0);
    // Short-Wave is stored in QrSolar. 
    if(this->db().foundObject<volScalarField>("qrSolar"))
      {      
	qSW = patch.lookupPatchField<volScalarField, scalar>("qrSolar");
      }    
    // Long-Wave (Built-Environment)
    scalarField qLW(this->size(), 0.0);
    if(this->db().foundObject<volScalarField>("qr"))
      {      
	qLW = patch.lookupPatchField<volScalarField, scalar>("qr");
      }
    if(this->db().foundObject<volScalarField>("Qr"))
      {      
	qLW = patch.lookupPatchField<volScalarField, scalar>("Qr");
      }
    // Long-Wave (sky)
    vector verticalDir(0,0,1);
    scalarField cosEpsilon=(verticalDir & -this->patch().nf());
    // Incoming Radiation from Clouds (not accounted in the exchanges)
    //    scalarField Tw_=patch.lookupPatchField<volScalarField, scalar>("T");
    scalarField qSky_=0.5*(1+cosEpsilon)*(213+5.5*(this->patchInternalField()-273.15));
    return (qSW+qLW+qSky_);
  }
  
  tmp<scalarField> ABLMOWallTFvPatchScalarField::latentHeat(scalarField& ustar,scalarField& nhat_,scalar MOL) const 
  {
    volScalarField qvTemp=this->db().lookupObject<volScalarField>("qv");
    const label patchi = this->patch().index();
    scalar kappa_=0.41;
    scalarField qvCell=qvTemp.boundaryField()[patchi].patchInternalField();    
    scalarField qvFace=this->patch().lookupPatchField<volScalarField, scalar>("qv");
    scalarField Tw_=this->patch().lookupPatchField<volScalarField, scalar>("T");
    scalar mavail=0.1;
    scalarField faceesat=610.94*exp(17.625*(Tw_-273.15)/(Tw_-30.11));
    scalarField faceqsat=(0.622 *faceesat/(101325-0.378*faceesat));
    // During spin-up latent heat can become too much
    // WRF suggests using a limit to avoid really high coefficient of moisture 
    scalarField denomq=max((log((nhat_+z0_)/z0_)-psiH(nhat_,MOL)+SMALL),2.0+0*Tw_);
    scalarField qstar=(faceqsat-qvCell)*kappa_/denomq;
    scalarField LH=1.225*2.2e6*mavail*ustar*qstar;
    minMax(LH);
    return (pos(LH)*LH);	  
  }
  
  tmp<scalarField> ABLMOWallTFvPatchScalarField::groundHeat() const
  {
    scalarField Tw=this->patch().lookupPatchField<volScalarField, scalar>("T");
    return (0.5*thermalCon_[0]/thickness_[0]*(Tw - T1_));
  }
  tmp<scalarField> ABLMOWallTFvPatchScalarField::firstLayer(scalarField& Twall,scalar dt_) const
  {
    scalar kWest=thermalCon_[0];
    scalar deltaxWest=0.5*thickness_[0];
    scalar kEast=0.5*(thermalCon_[0]+thermalCon_[1]);
    scalar deltaxEast=0.5*(thickness_[0]+thickness_[1]);
    scalar Coeff_=matDensity_[0]*specificHeat_[0]/dt_+kWest/deltaxWest+kEast/deltaxEast;
    scalarField Tw=this->patch().lookupPatchField<volScalarField, scalar>("T");
    return (1/Coeff_*(kWest/deltaxWest*Tw+matDensity_[0]*specificHeat_[0]/dt_*T1_+kEast/deltaxEast*T2_));
  }

  tmp<scalarField> ABLMOWallTFvPatchScalarField::secondLayer(scalar dt_) const
  {
    scalar kWest=0.5*(thermalCon_[0]+thermalCon_[1]);
    scalar deltaxWest=0.5*(thickness_[0]+thickness_[1]);
    scalar kEast=0.5*(thermalCon_[1]+thermalCon_[2]);
    scalar deltaxEast=0.5*(thickness_[1]+thickness_[2]);
    scalar Coeff_=matDensity_[1]*specificHeat_[1]/dt_+kWest/deltaxWest+kEast/deltaxEast;
    return (1/Coeff_*(kWest/deltaxWest*T1_+matDensity_[1]*specificHeat_[1]/dt_*T2_+kEast/deltaxEast*T3_));   
  }

  tmp<scalarField> ABLMOWallTFvPatchScalarField::thirdLayer(scalar dt_) const
  {
    scalar kWest=0.5*(thermalCon_[1]+thermalCon_[2]);
    scalar deltaxWest=0.5*(thickness_[1]+thickness_[2]);
    scalar kEast=0.5*(thermalCon_[2]+thermalCon_[3]);
    scalar deltaxEast=0.5*(thickness_[2]+thickness_[3]);
    scalar Coeff_=matDensity_[2]*specificHeat_[2]/dt_+kWest/deltaxWest+kEast/deltaxEast;
    return (1/Coeff_*(kWest/deltaxWest*T2_+matDensity_[2]*specificHeat_[2]/dt_*T3_+kEast/deltaxEast*T4_));
  }

  tmp<scalarField> ABLMOWallTFvPatchScalarField::fourthLayer(scalar dt_) const
  {
    scalar kWest=0.5*(thermalCon_[2]+thermalCon_[3]);
    scalar deltaxWest=0.5*(thickness_[2]+thickness_[3]);
    scalar kEast=thermalCon_[3];
    scalar deltaxEast=0.5*thickness_[3];
    scalar Coeff_=matDensity_[3]*specificHeat_[3]/dt_+kWest/deltaxWest+kEast/deltaxEast;
    return (1/Coeff_*(kWest/deltaxWest*T3_+matDensity_[3]*specificHeat_[3]/dt_*T4_+kEast/deltaxEast*Tdeep_));   
  }
  void ABLMOWallTFvPatchScalarField::minMax(scalarField& field) const
  {
    scalar qMin_=min(field);                                                                                                
    scalar qMax_=max(field);                                                                                                
    reduce(qMin_, minOp<scalar>());                                                                                         
    reduce(qMax_, maxOp<scalar>());                                                                                         
    //    Info<<this->patch().name()<<" Field  "<<qMin_<<"  "<<qMax_<<endl;       
  }
  void ABLMOWallTFvPatchScalarField::write(Ostream& os) const
  {
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("z0")<< z0_ << token::END_STATEMENT << nl;
    matDensity_.writeEntry("matDensity",os);
    specificHeat_.writeEntry("specificHeat",os);
    thermalCon_.writeEntry("thermalCon", os);
    thickness_.writeEntry("thickness", os);
    T1_.writeEntry("T1", os);
    T2_.writeEntry("T2", os);
    T3_.writeEntry("T3", os);
    T4_.writeEntry("T4", os);
    Tdeep_.writeEntry("Tdeep", os);   
    //    writeEntry("value", os);

  }
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  makePatchTypeField
  (
   fvPatchScalarField,
   ABLMOWallTFvPatchScalarField
   );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// ************************************************************************* //
