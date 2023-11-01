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

#include "ABLLouisWallTFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"
#include "scalarIOList.H"
#include "boundaryRadiationProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


 

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  ABLLouisWallTFvPatchScalarField::
  ABLLouisWallTFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(p, iF),
    z0_(0.1),
    mavail_(0.1),
    matDensity_(),
    specificHeat_(),
    thermalCon_(),
    thickness_()
  {}


  ABLLouisWallTFvPatchScalarField::
  ABLLouisWallTFvPatchScalarField
  (
   const ABLLouisWallTFvPatchScalarField& ptf,
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const fvPatchFieldMapper& mapper
   )
    :
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    z0_(ptf.z0_),
    mavail_(ptf.mavail_),
    matDensity_(ptf.matDensity_),
    specificHeat_(ptf.specificHeat_),
    thermalCon_(ptf.thermalCon_),
    thickness_(ptf.thickness_),
    sensibleHeat_(),
    Tw_(),
    T1_(),
    T2_(),
    T3_(),
    T4_(),
    Tdeep_()
  {
    sensibleHeat_.setSize(mapper.size());
    sensibleHeat_.map(ptf.sensibleHeat_,mapper);
    Tw_.setSize(mapper.size());
    Tw_.map(ptf.Tw_, mapper);
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
  ABLLouisWallTFvPatchScalarField::
  ABLLouisWallTFvPatchScalarField
  (
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const dictionary& dict
   )
    :
    fixedValueFvPatchScalarField(p, iF),
    z0_(dict.lookupOrDefault<scalar>("z0",0.1)),
    mavail_(dict.lookupOrDefault<scalar>("mavail",0.1)),
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
    if(dict.found("sensibleHeat"))
      {
	sensibleHeat_=scalarField("sensibleHeat",dict,p.size());
      }
    if (dict.found("Tw"))
      {
        Tw_= scalarField("Tw", dict, p.size());
      }   
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


  ABLLouisWallTFvPatchScalarField::
  ABLLouisWallTFvPatchScalarField
  (
   const ABLLouisWallTFvPatchScalarField& rwfpsf
   )
    :
    fixedValueFvPatchScalarField(rwfpsf),
    z0_(rwfpsf.z0_),
    mavail_(rwfpsf.mavail_),
    matDensity_(rwfpsf.matDensity_),
    specificHeat_(rwfpsf.specificHeat_),
    thermalCon_(rwfpsf.thermalCon_),
    thickness_(rwfpsf.thickness_),
    sensibleHeat_(rwfpsf.sensibleHeat_),
    Tw_(rwfpsf.Tw_),    
    T1_(rwfpsf.T1_),
    T2_(rwfpsf.T2_),
    T3_(rwfpsf.T3_),
    T4_(rwfpsf.T4_),
    Tdeep_(rwfpsf.Tdeep_)   
  {}



  ABLLouisWallTFvPatchScalarField::
  ABLLouisWallTFvPatchScalarField
  (
   const ABLLouisWallTFvPatchScalarField& rwfpsf,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    fixedValueFvPatchScalarField(rwfpsf, iF),
    z0_(rwfpsf.z0_),
    mavail_(rwfpsf.mavail_),
    matDensity_(rwfpsf.matDensity_),
    specificHeat_(rwfpsf.specificHeat_),    
    thermalCon_(rwfpsf.thermalCon_),
    thickness_(rwfpsf.thickness_),
    sensibleHeat_(rwfpsf.sensibleHeat_),
    Tw_(rwfpsf.Tw_),    
    T1_(rwfpsf.T1_),
    T2_(rwfpsf.T2_),
    T3_(rwfpsf.T3_),
    T4_(rwfpsf.T4_),
    Tdeep_(rwfpsf.Tdeep_)    

  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void ABLLouisWallTFvPatchScalarField::autoMap
  (
   const fvPatchFieldMapper& m
   )
  {
    fixedValueFvPatchScalarField::autoMap(m);
    sensibleHeat_.autoMap(m);    
    Tw_.autoMap(m);    
    T1_.autoMap(m);
    T2_.autoMap(m);
    T3_.autoMap(m);
    T4_.autoMap(m);
    Tdeep_.autoMap(m);
  }


  void ABLLouisWallTFvPatchScalarField::rmap
  (
   const fvPatchScalarField& ptf,
   const labelList& addr
   )
  {
    fixedValueFvPatchScalarField::rmap(ptf, addr);
    const ABLLouisWallTFvPatchScalarField& ewhftpsf =
      refCast<const ABLLouisWallTFvPatchScalarField>(ptf);
    sensibleHeat_.rmap(ewhftpsf.sensibleHeat_, addr);    
    Tw_.rmap(ewhftpsf.Tw_, addr);
    T1_.rmap(ewhftpsf.T1_, addr);
    T2_.rmap(ewhftpsf.T2_, addr);
    T3_.rmap(ewhftpsf.T3_, addr);
    T4_.rmap(ewhftpsf.T4_, addr);
    Tdeep_.rmap(ewhftpsf.Tdeep_, addr);
    return;
  }

  // * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

  void ABLLouisWallTFvPatchScalarField::updateCoeffs() 
  { 
    if (this->updated())
      {
	return;
      }
    // Initial Start up 
    scalarField field=patch().lookupPatchField<volScalarField, scalar>("qr"); 
    scalar qMin_=min(field);
    scalar qMax_=max(field);
    reduce(qMin_, minOp<scalar>());
    reduce(qMax_, maxOp<scalar>());
    if(qMin_==qMax_)
      {
	scalarField::operator=(Tw_);
	fixedValueFvPatchScalarField::updateCoeffs();
	return;
      }
    // Read the driving variables from the solver 
    const scalarIOList& Tvalues=this->db().lookupObject<scalarIOList>("timeProfile");
    const scalarIOList& MOvalues=this->db().lookupObject<scalarIOList>("MOProfile");
    const scalarIOList& TRefvalues=this->db().lookupObject<scalarIOList>("TProfile");
    const scalarIOList& PBLHvalues=this->db().lookupObject<scalarIOList>("PBLHProfile");
    scalarField timevalues(Tvalues.size(),0.0);
    scalarField trefvalues(Tvalues.size(),0.0);
    scalarField movalues(Tvalues.size(),0.0);
    scalarField pblhvalues(Tvalues.size(),0.0);
    forAll(Tvalues,i)
      {
        timevalues[i] = Tvalues[i];
        trefvalues[i]=TRefvalues[i];
	movalues[i] = MOvalues[i];
        pblhvalues[i] = PBLHvalues[i];
      }
    const fvPatch& patch=this->patch();
    const label patchi = patch.index();     
    const volVectorField& U=this->db().lookupObject<volVectorField>("U");
    const vectorField  Uc = U.boundaryField()[patchi].patchInternalField();
    scalar kappa_=0.41;
    scalar MOL=interpolateXY(U.time().value(),timevalues,movalues);
    if(mag(MOL)>25)
      MOL=MOL;
    else if(neg(MOL) && MOL> -25)
      MOL=-25;
    else if(pos(MOL) && mag(MOL)<25)
      MOL=25;
    scalar PBLH=interpolateXY(U.time().value(),timevalues,pblhvalues);
    PBLH=1000;
    scalar TRef=interpolateXY(U.time().value(),timevalues,trefvalues);
    // Adjust wall distance 
    vectorField normalV_=patch.nf();
    scalarField y_=mag((patch.Cn()-patch.Cf()) & patch.nf());
    scalarField nhat_=max(y_,0.01+0*y_)+z0_;  //Avoid thin gap blow up 
    // Extract the parallel velocity 
    vectorField UpParallel=Uc - (Uc & normalV_)*normalV_;
    // Radiation 
    scalarField radiation_=patch.lookupPatchField<volScalarField, scalar>("qr");
    minMax("Radiation:",radiation_);
    scalar overflow=0;
    scalar underflow=0;
    scalar totalface=0;
    // Shield the ground with trees
    const volScalarField treeTemp=this->db().lookupObject<volScalarField>("treeBlanking");
    const scalarField treeCell=treeTemp.boundaryField()[patchi].patchInternalField();
    radiation_=(1-pos(treeCell))*radiation_;
    // convective velocity wstar from values at last step for unstable boundary layer 
    sensibleHeat_=min(sensibleHeat_,0*sensibleHeat_);
    scalarField wstar=pow(-9.81/this->patchInternalField()*PBLH*neg(MOL)*sensibleHeat_,0.33);
    // Louis Formula for sensible heat flux
    scalarField deltaTheta=this->patchInternalField()-Tw_;
    scalarField louisA=kappa_/log((nhat_+z0_)/z0_);  
    scalar Cstar=7.4;
    scalarField magUp=mag(UpParallel);
    magUp=max(magUp,0.1+0*magUp);
    scalarField z01_=z0_+0*magUp;
    scalarField ustar=louisA*magUp*sqrt(louisF(deltaTheta,this->patchInternalField(),magUp,nhat_,MOL,Cstar,z01_)); 
    Cstar=5.3;
    scalar Czil_=pow(10,-kappa_*z0_/0.07);
    scalar nu=1e-5;
    scalarField zt_=z0_*exp(-kappa_*Czil_*sqrt(ustar*z0_/nu));
    louisA=kappa_/log((nhat_+zt_)/zt_); 
    scalarField ustarThetaStar=louisA/0.7*deltaTheta*sqrt(louisF(deltaTheta,this->patchInternalField(),magUp,nhat_,MOL,Cstar,zt_))*sqrt(pow(ustar,2)+pow(wstar,2));
    sensibleHeat_=1225*ustarThetaStar;
    minMax("u*:",ustar);
    minMax("SH:",sensibleHeat_);
    // Latent Heat from Dudhia
    const volScalarField qvTemp=this->db().lookupObject<volScalarField>("qv");
    const scalarField qVCell=qvTemp.boundaryField()[patchi].patchInternalField();
    scalar ka=2.4e-5;
    scalarField den=log(kappa_*ustar*nhat_/ka+nhat_/0.01);
    scalarField faceesat=610.94*exp(17.625*(Tw_-273.15)/(Tw_-30.11));
    scalarField qvFace=(0.622 *faceesat/(101325-0.378*faceesat));
    scalarField qstar=(pos(treeCell)*0.5+neg(treeCell)*mavail_)*kappa_*(qVCell-qvFace)/den;
    scalarField latentHeat_=1.225*2.2e6*min(ustar,0.5+0*ustar)*qstar;
    minMax("LH:",latentHeat_);
    // SEB Update 
    //    if((word(U.mesh().ddtScheme("default")) != "steadyState"))                     
    //      dt_=U.mesh().time().deltaTValue();
    scalarField Twall=this->patchInternalField();
    scalarField groundHeat_=(groundHeat());
    //minMax("G:",groundHeat_);
    // Correct for Source Term from solver 
    const volScalarField sourceBlanking=this->db().lookupObject<volScalarField>("sourceBlanking");
    const scalarField sourceBlankingValue=sourceBlanking.boundaryField()[patchi].patchInternalField();
    const volScalarField temperatureSource=this->db().lookupObject<volScalarField>("temperatureSource");
    const scalarField temperatureSourceValue=temperatureSource.boundaryField()[patchi].patchInternalField();
    scalarField externalSource=0*groundHeat_;
    // If there is a source term no ground-heating can happen 
    forAll(externalSource,faceI)
      {
	if(sourceBlankingValue[faceI]>0)
	  {
	    externalSource[faceI]=temperatureSourceValue[faceI];
	    groundHeat_[faceI]=0.0;
	  }
      }
    // Temporary place-holder for water body radiation penetration 
    scalar radiationExtinct_=1.0;
    if(mavail_>0.9)
      radiationExtinct_=0.4;
    scalarField residual=(radiationExtinct_*radiation_+sensibleHeat_+latentHeat_+groundHeat_+externalSource);
    //minMax("R:",residual);
    scalar dt_=600;
    scalar implicitTerm=matDensity_[0]*specificHeat_[0]/dt_; 
    Twall=1.0/implicitTerm*(Tw_*implicitTerm+residual);
    scalarField Tc=this->patchInternalField();
    /*forAll(Tw_,faceI)
      {
	if(Tw_[faceI]<297)
	  Pout<<patch.name()<<"T:"<<Tw_[faceI]<<" Rad:"<<radiation_[faceI]
	      <<" SH:"<<sensibleHeat_[faceI]<<" "<<ustar[faceI]<<
	    " LH:"<<latentHeat_[faceI]<<" G:"<<groundHeat_[faceI]<<" "<<radiationExtinct_<<endl;
	    }*/ 
    // Correct for Radiation going crazy
    overflow=0; 
    underflow=0;
    if((word(U.mesh().ddtScheme("default")) == "steadyState"))
      {
	forAll(radiation_,faceI)
	  {
	    if(Tw_[faceI]<(TRef-3) ||Tw_[faceI]>(TRef+8))
	      {
		scalar T_=0;
		vector faceCor=patch.Cf()[faceI];
		scalar counter=0;
		forAll(Twall,faceJ)
		  {
		    scalar distance=mag(faceCor-patch.Cf()[faceJ]);
		    if(distance>1e-5 && distance < 2 )
		      {
			T_=T_+Tw_[faceJ];
			counter=counter+1;
		      }
		  }
		if(counter>0)
		  {
		    Twall[faceI]=T_/counter;
		  }
		else
		  {
		    Twall[faceI]=Tc[faceI];
		  }
		if(Tw_[faceI]<(TRef-3))
		  underflow=underflow+1;
		if(Tw_[faceI]>(TRef+8))
		  overflow=overflow+1;
	      }
	  }
	reduce(overflow,sumOp<scalar>());                                                         
	reduce(underflow,sumOp<scalar>());                                                       
	if(overflow>0)
	  Info<<"Cooling:"<<patch.name()<<"  "<<overflow<<" out of "<<totalface<<endl;
	if(underflow>0)
	  Info<<"Heating:"<<patch.name()<<"  "<<underflow<<" out of "<<totalface<<endl;
      }
    // Store the virtual heat in sensible heat for next step 	
    // Place holder for future immersed boundary 
    /*if(this->db().foundObject<volScalarField>("buildingBlanking"))
      {
      const volScalarField building=this->db().lookupObject<volScalarField>("buildingBlanking");
      const scalarField buildingBlanking=building.boundaryField()[patchi].patchInternalField();
      Twall=(1-buildingBlanking)*this->patchInternalField()+buildingBlanking*Twall;
      } */
    // Avoid the residual calculation error from slowing down convergence 
    scalar cutOff=0.008;
    scalarField zone=0*Tw_;
    underflow=0;
    forAll(Twall,faceI)
      {
	if((mag(Twall[faceI]-Tw_[faceI])<cutOff)||zone[faceI]==1)
	  {
	    Twall[faceI]=Tw_[faceI];
	    underflow=underflow+1;
	  }
      }
    sensibleHeat_=(1+0.61*qVCell)*ustarThetaStar+0.61*this->patchInternalField()*ustar*qstar;
    //reduce(underflow,sumOp<scalar>());    
    //if(underflow>0)
    //  Info<<"Residual correction:"<<patch.name()<<"  "<<underflow<<" out of "<<totalface<<endl;
                                                     
    // Relaxed Update to match SIMPLE 
    Tw_=0.5*Tw_+0.5*Twall;    
    //minMax("DT:",mag(Twall-Tw_));
    minMax("Wall Temperature:",Tw_);
    //minMax("Tc:",Tc);
    // Update Temperature
    // Wall heat conduction
    //    if((word(U.mesh().ddtScheme("default")) != "steadyState"))
    //      dt_=U.mesh().time().deltaTValue();
    scalar materialDiffusion=thermalCon_[0]/(matDensity_[0]*specificHeat_[0]);
    scalar dtConduction=0.5*thickness_[0]*thickness_[0]/materialDiffusion;
    scalar steps=dt_/dtConduction;
    scalar finalsteps=0;
    //Info<<patch.name()<<"Conduction time step:"<<dtConduction<<endl;
    if(dtConduction>dt_)
      {
	T1_=firstLayer(dt_);
	// Second-layer
	T2_=secondLayer(dt_);
	// Third-layer
	T3_=thirdLayer(dt_);
	// Last-Layer
	T4_=fourthLayer(dt_);
      }
    else
      {
	//Info<<patch.name()<<"Conduction steps:"<<steps<<endl;
	while(finalsteps<steps)
	  {
	    //Info<<patch.name()<<finalsteps<<endl;
	    T1_=firstLayer(dtConduction);
	    // Second-layer
	    T2_=secondLayer(dtConduction);
	    // Third-layer
	    T3_=thirdLayer(dtConduction);
	    // Last-Layer
	    T4_=fourthLayer(dtConduction);
	    finalsteps=finalsteps+1;
	  }
      }
    scalarField::operator=(Tw_);
    fixedValueFvPatchScalarField::updateCoeffs();
  }

  tmp<scalarField> ABLLouisWallTFvPatchScalarField::louisF(scalarField& deltaTheta,const scalarField& Tcell,scalarField& Uc,const scalarField& nhat_,scalar& MOL, scalar& louisCstar,scalarField& zt_) const
  {
    // Replace zt_ with zt_
    // Unstable
    scalarField Rib=9.81*nhat_*deltaTheta/(Tcell*sqr(Uc)+SMALL);
    Rib=max(Rib,-10+0*Rib);
    scalar louisB=9.4;
    scalar kappa_=0.41;
    scalarField louisC=louisCstar*louisB*sqr(kappa_/log((nhat_+zt_)/zt_))*sqrt(nhat_/zt_);
    scalar louisBstar=4.7;
    return neg(Rib)*(1-louisB*Rib/(1+louisC*sqrt(mag(Rib))))+pos(Rib)/sqr(1+louisBstar*Rib);
  }

  
  tmp<scalarField> ABLLouisWallTFvPatchScalarField::groundHeat() const
  {
    return (0.5*thermalCon_[0]/thickness_[0]*(T1_-Tw_));
  }
  tmp<scalarField> ABLLouisWallTFvPatchScalarField::firstLayer(scalar dt_) const
  {
    scalar kWest=thermalCon_[0];
    scalar deltaxWest=0.5*thickness_[0];
    scalar kEast=0.5*(thermalCon_[0]+thermalCon_[1]);
    scalar deltaxEast=0.5*(thickness_[0]+thickness_[1]);
    scalar Coeff_=matDensity_[0]*specificHeat_[0]/dt_+kWest/deltaxWest+kEast/deltaxEast;
    return (1/Coeff_*(kWest/deltaxWest*Tw_+matDensity_[0]*specificHeat_[0]/dt_*T1_+kEast/deltaxEast*T2_));
  }

  tmp<scalarField> ABLLouisWallTFvPatchScalarField::secondLayer(scalar dt_) const
  {
    scalar kWest=0.5*(thermalCon_[0]+thermalCon_[1]);
    scalar deltaxWest=0.5*(thickness_[0]+thickness_[1]);
    scalar kEast=0.5*(thermalCon_[1]+thermalCon_[2]);
    scalar deltaxEast=0.5*(thickness_[1]+thickness_[2]);
    scalar Coeff_=matDensity_[1]*specificHeat_[1]/dt_+kWest/deltaxWest+kEast/deltaxEast;
    return (1/Coeff_*(kWest/deltaxWest*T1_+matDensity_[1]*specificHeat_[1]/dt_*T2_+kEast/deltaxEast*T3_));   
  }

  tmp<scalarField> ABLLouisWallTFvPatchScalarField::thirdLayer(scalar dt_) const
  {
    scalar kWest=0.5*(thermalCon_[1]+thermalCon_[2]);
    scalar deltaxWest=0.5*(thickness_[1]+thickness_[2]);
    scalar kEast=0.5*(thermalCon_[2]+thermalCon_[3]);
    scalar deltaxEast=0.5*(thickness_[2]+thickness_[3]);
    scalar Coeff_=matDensity_[2]*specificHeat_[2]/dt_+kWest/deltaxWest+kEast/deltaxEast;
    return (1/Coeff_*(kWest/deltaxWest*T2_+matDensity_[2]*specificHeat_[2]/dt_*T3_+kEast/deltaxEast*T4_));
  }

  tmp<scalarField> ABLLouisWallTFvPatchScalarField::fourthLayer(scalar dt_) const
  {
    scalar kWest=0.5*(thermalCon_[2]+thermalCon_[3]);
    scalar deltaxWest=0.5*(thickness_[2]+thickness_[3]);
    scalar kEast=thermalCon_[3];
    scalar deltaxEast=0.5*thickness_[3];
    scalar Coeff_=matDensity_[3]*specificHeat_[3]/dt_+kWest/deltaxWest+kEast/deltaxEast;
    return (1/Coeff_*(kWest/deltaxWest*T3_+matDensity_[3]*specificHeat_[3]/dt_*T4_+kEast/deltaxEast*Tdeep_));   
  }
  void ABLLouisWallTFvPatchScalarField::minMax(word variableName,scalarField& field) const
  {
    scalar qMin_=min(field);                                                                                        
    scalar qMax_=max(field);
    reduce(qMin_, minOp<scalar>());     
    reduce(qMax_, maxOp<scalar>());
    Info<<this->patch().name()<<" "<<variableName<<":[ "<<qMin_<<"  "<<qMax_<<" ]"<<endl;       
  }
  void ABLLouisWallTFvPatchScalarField::minMax(word variableName,const scalarField& field) const
  {
    scalar qMin_=min(field);
    scalar qMax_=max(field);
    reduce(qMin_, minOp<scalar>());
    reduce(qMax_, maxOp<scalar>());                                                                 
    Info<<this->patch().name()<<" "<<variableName<<":[ "<<qMin_<<"  "<<qMax_<<" ]"<<endl;       
  }

  void ABLLouisWallTFvPatchScalarField::write(Ostream& os) const
  {
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("z0")<< z0_ << token::END_STATEMENT << nl;
    os.writeKeyword("mavail")<< mavail_ << token::END_STATEMENT << nl;
    matDensity_.writeEntry("matDensity",os);
    specificHeat_.writeEntry("specificHeat",os);
    thermalCon_.writeEntry("thermalCon", os);
    thickness_.writeEntry("thickness", os);
    sensibleHeat_.writeEntry("sensibleHeat",os);
    Tw_.writeEntry("Tw", os);
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
   ABLLouisWallTFvPatchScalarField
   );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

  // ************************************************************************* //
