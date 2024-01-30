/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | www.openfoam.com
  \\/     M anipulation  |
  -------------------------------------------------------------------------------
  Copyright (C) 2011-2017 OpenFOAM Foundation
  Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "wideBandDiffusiveRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOM.H"
#include "wideBandAbsorptionEmission.H"
#include "constants.H"
#include "boundaryRadiationProperties.H"
#include "wallFvPatch.H"
// Added for Sky
#include "interpolateXY.H"
#include "scalarIOList.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF
 )
  :
  mixedFvPatchScalarField(p, iF)
{
  refValue() = 0.0;
  refGrad() = 0.0;
  valueFraction() = 1.0;
}


Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
 const wideBandDiffusiveRadiationMixedFvPatchScalarField& ptf,
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF,
 const fvPatchFieldMapper& mapper
 )
  :
  mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF,
 const dictionary& dict
 )
  :
  mixedFvPatchScalarField(p, iF)
{
  if (dict.found("refValue"))
    {
      fvPatchScalarField::operator=
        (
	 scalarField("value", dict, p.size())
	 );
      refValue() = scalarField("refValue", dict, p.size());
      refGrad() = scalarField("refGradient", dict, p.size());
      valueFraction() = scalarField("valueFraction", dict, p.size());
    }
  else
    {
      refValue() = 0.0;
      refGrad() = 0.0;
      valueFraction() = 1.0;

      fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
 const wideBandDiffusiveRadiationMixedFvPatchScalarField& ptf
 )
  :
  mixedFvPatchScalarField(ptf)
{}


Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
 const wideBandDiffusiveRadiationMixedFvPatchScalarField& ptf,
 const DimensionedField<scalar, volMesh>& iF
 )
  :
  mixedFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
updateCoeffs()
{
  if (this->updated())
    {
      return;
    }

  // Since we're inside initEvaluate/evaluate there might be processor
  // comms underway. Change the tag we use.
  int oldTag = UPstream::msgType();
  UPstream::msgType() = oldTag+1;
  const scalarField& Tp=patch().lookupPatchField<volScalarField, scalar>("T");
  const radiationModel& radiation =
    db().lookupObject<radiationModel>("radiationProperties");

  const fvDOM& dom(refCast<const fvDOM>(radiation));

  label rayId = -1;
  label lambdaId = -1;
  dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

  const label patchi = patch().index();

  if (dom.nLambda() == 0)
    {
      FatalErrorInFunction
	<< " a non-grey boundary condition is used with a grey "
	<< "absorption model" << nl << exit(FatalError);
    }

  scalarField& Iw = *this;
  const vectorField n(patch().Sf()/patch().magSf());

  radiativeIntensityRay& ray =
    const_cast<radiativeIntensityRay&>(dom.IRay(rayId));

  const scalarField nAve(n & ray.dAve());

  //  ray.qr().boundaryFieldRef()[patchi] += Iw*nAve;

  const scalarField Eb
    (
     dom.blackBody().bLambda(lambdaId).boundaryField()[patchi]
     );

  const boundaryRadiationProperties& boundaryRadiation =
    boundaryRadiationProperties::New(internalField().mesh());


  const tmp<scalarField> temissivity
    (
     boundaryRadiation.emissivity(patch().index(), lambdaId)
     );
  const scalarField& emissivity = temissivity();
  const tmp<scalarField> tabsorptivity
    (
     boundaryRadiation.absorptivity(patch().index(), lambdaId)
     );
  const scalarField& absorptivity = tabsorptivity();
  //  Pout<<patch().name()<<max(absorptivity)<<endl;
  const tmp<scalarField> ttransmissivity
    (
     boundaryRadiation.transmissivity(patch().index(), lambdaId)
     );
  const scalarField& transmissivity = ttransmissivity();
  scalarField& qem = ray.qem().boundaryFieldRef()[patchi];
  scalarField& qin = ray.qin().boundaryFieldRef()[patchi];
  scalarField& qcb = ray.qcb().boundaryFieldRef()[patchi];
  // Looking for the ray closest to the Sun direction                                                                         
  label SunRayId(-1);
  scalar maxSunRay = -GREAT;
  const vector sunDir = dom.solarCalc().direction();
  for (label rayI=0; rayI < dom.nRay(); rayI++)
    {
      const vector& iD = dom.IRay(rayI).d();
      scalar dir = sunDir & iD;
      if (dir > maxSunRay)
        {
          maxSunRay = dir;
          SunRayId = rayI;
        }
    }
  scalarField blocking=1+0*absorptivity;
  if(rayId == SunRayId)
    {
      blocking=absorptivity;
      ray.qr().boundaryFieldRef()[patchi] = blocking*Iw*nAve;
    }
  else
    {
      scalarField tempQr=Iw*nAve;
      /*scalar Iwcounter=0;
      forAll(tempQr,faceI)
	{
	  if(tempQr[faceI]< -400)
	    {
	      tempQr[faceI]=-400.0;
	      Iwcounter=Iwcounter+1;
	    }
	  if(tempQr[faceI]> 800)
            {
              tempQr[faceI]=800.0;
              Iwcounter=Iwcounter+1;
            }

	}
      reduce(Iwcounter,sumOp<scalar>());
      if(Iwcounter>0)
	{
	  Info<<"DOM diverged for ray:"<<rayId<<" on patch:"<<this->patch().name()
	      <<" at "<<Iwcounter<<" points"<<endl;
	}
      /*scalar qMin_=min(tempQr);
      scalar qMax_=max(tempQr);
      reduce(qMin_, minOp<scalar>());
      reduce(qMax_, maxOp<scalar>());
      if(qMin_< -200 || qMax_ > 800)
	{
	  Info<<"DOM diverged for ray:"<<rayId<<" on patch:"<<this->patch().name()
	      <<" and neglecting contribution"<<"   "<<qMin_<<"   "<<qMax_<<endl;
	  ray.qr().boundaryFieldRef()[patchi] = 0*Iw*nAve;
	}
	else */
      //ray.qr().boundaryFieldRef()[patchi] = Iw*nAve;
      ray.qr().boundaryFieldRef()[patchi] = tempQr;
    }
  scalar qMin_=min(blocking*Iw*nAve);			\
  scalar qMax_=max(blocking*Iw*nAve);
  reduce(qMin_, minOp<scalar>());
  reduce(qMin_, minOp<scalar>());
  if(qMin_< - 500 || qMax_> 1000)
    {
     
      Info<<"Patch:"<<patch().name();
      Info<<" Ray:"<<rayId<<"  "<<SunRayId<<"  ";
      minMax("Ray:", blocking*Iw*nAve);
    }
  // Calculate Ir into the wall on the same lambdaId
  scalarField Ir(patch().size(), Zero);
  forAll(Iw, facei)
    {
      for (label rayi=0; rayi < dom.nRay(); rayi++)
        {
	  const vector& d = dom.IRay(rayi).d();

	  if ((-n[facei] & d) < 0.0)
            {
	      // q into the wall
	      const scalarField& IFace =
		dom.IRay(rayi).ILambda(lambdaId).boundaryField()[patchi];

	      const vector& rayDave = dom.IRay(rayi).dAve();
	      Ir[facei] += IFace[facei]*(n[facei] & rayDave);
            }
        }
    }
  scalarField qSky_=0*Ir;
  //Info<<"Sky Radiation:"<<cSky_<<endl;
  if (dom.useSolarLoad())
    {
      const volScalarField& qSolar =this->db().lookupObject<volScalarField>("qrSolar");
      Ir += qSolar.boundaryField()[patch().index()];
      if (isA<wallFvPatch>(patch()))
	{
	  vector verticalDir(0,0,1);
	  scalarField cosEpsilon=(verticalDir & -this->patch().nf());
	  const volVectorField& U=this->db().lookupObject<volVectorField>("U");
	  const scalarIOList& Tvalues=this->db().lookupObject<scalarIOList>("timeProfile");
	  const scalarIOList& TRefvalues=this->db().lookupObject<scalarIOList>("TProfile");
	  scalarField timevalues(Tvalues.size(),0.0);
	  scalarField trefvalues(Tvalues.size(),0.0);
	  forAll(Tvalues,i)
	    {
	      timevalues[i] = Tvalues[i];
	      trefvalues[i]=TRefvalues[i];
	    }
	    scalar TRef=interpolateXY(U.time().value(),timevalues,trefvalues); 
	  scalar Ta=TRef-13;
	  qSky_=0.5*(1+cosEpsilon)*emissivity*5.67e-8*(pow(Ta,4)-pow(Tp,4));
	  Ir+=qSky_;
	}
    }

  scalarField Iexternal(this->size(), 0.0);
  if (dom.useExternalBeam())
    {
      const vector sunDir = dom.solarCalc().direction();
      const scalar directSolarRad =
	dom.solarCalc().directSolarRad()
	*dom.spectralDistribution()[lambdaId];

      //label nRaysBeam = dom.nRaysBeam();
      //label SunRayId(-1);
      //      scalar maxSunRay = -GREAT;

      // Looking for the ray closest to the Sun direction
      for (label rayI=0; rayI < dom.nRay(); rayI++)
        {
	  const vector& iD = dom.IRay(rayI).d();
	  scalar dir = sunDir & iD;
	  if (dir > maxSunRay)
            {
	      maxSunRay = dir;
	      SunRayId = rayI;
            }
        }
      
      if (rayId == SunRayId)
        {
	  const scalarField nAve(n & dom.IRay(rayId).dAve());
	  forAll(Iexternal, faceI)
            {
	      Iexternal[faceI] = directSolarRad/mag(dom.IRay(rayId).dAve());
            }  
	  if (isA<wallFvPatch>(patch()))
	    {
	      forAll(Iw, facei)
		{
		  const vector& d = dom.IRay(rayId).d();
		  if ((-n[facei] & d) > 0.0)
		    qcb[facei]=0.0;
		  else
		    qcb[facei]= absorptivity[facei]*Iw[facei]*nAve[facei];
		}
	    }
	  //minMax("SW:",qcb);
	}
    }


  if (isA<wallFvPatch>(patch()))
    {
      /*      if(patch().name()=="set64_28")
	{
	  Info<<"ID"<<rayId<<"  "<<SunRayId<<endl;
	  minMax("qr:",blocking*Iw*nAve);
	  minMax("Incoming:",Iw*nAve);
	  minMax("NonEb:",Ir*(1.0 - emissivity));
	  minMax("Eb:",emissivity*Eb);
	  }*/
      forAll(Iw, facei)
	{
	  const vector& d = dom.IRay(rayId).d();

	  if ((-n[facei] & d) > 0.0)
	    {
	      // direction out of the wall
	      refGrad()[facei] = 0.0;
	      valueFraction()[facei] = 1.0;
	      refValue()[facei] =
                Iexternal[facei]*transmissivity[facei]
		+ (
		   Ir[facei]*(1.0 - emissivity[facei])
		   + emissivity[facei]*physicoChemical::sigma.value()
		   * pow4(Tp[facei])
		   )/pi;

		   //		   + emissivity[facei]*Eb[facei]
		   //		   +qSky_[facei]
		   //		   )/pi;

	      // Emitted heat flux from this ray direction (sum over lambdaId)
	      qem[facei] += refValue()[facei]*nAve[facei];
	    }
	  else
	    {
	      // direction into the wall
	      valueFraction()[facei] = 0.0;
	      refGrad()[facei] = 0.0;
	      refValue()[facei] = 0.0; //not used

	      // Incident heat flux on this ray direction (sum over lambdaId)
	      qin[facei] += Iw[facei]*nAve[facei];
	    }
	}
    }
  else
    {
      // Adjust the emissivity for air 
      const volScalarField& T_=this->db().lookupObject<volScalarField>("T");
      const scalarField& Tw=T_.boundaryField()[patchi];     
      forAll(Iw, facei)
	{
	  const vector& d = dom.IRay(rayId).d();  
	  if ((-n[facei] & d) > 0.0)
	    {
	      // Swinbank's formula for emissivity 
	      scalar emissivity_air=9.365e-6*sqr(Tw[facei]);
	      // direction out of the wall
	      refGrad()[facei] = 0.0;
	      valueFraction()[facei] = 1.0;
	      refValue()[facei] =
                Iexternal[facei]*transmissivity[facei]
		+ (
		   Ir[facei]*(1.0 - emissivity[facei]*emissivity_air)
                   + emissivity[facei]*emissivity_air*physicoChemical::sigma.value()
                   * pow4(Tp[facei])
		   )/pi;

	      // Emitted heat flux from this ray direction (sum over lambdaId)
	      qem[facei] += refValue()[facei]*nAve[facei];
	    }
	  else
	    {
	      // direction into the wall
	      valueFraction()[facei] = 0.0;
	      refGrad()[facei] = 0.0;
	      refValue()[facei] = 0.0; //not used

	      // Incident heat flux on this ray direction (sum over lambdaId)
	      qin[facei] += Iw[facei]*nAve[facei];
	    }
	}
    }

  // Restore tag
  UPstream::msgType() = oldTag;

  mixedFvPatchScalarField::updateCoeffs();
}

  void Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::minMax(word variableName,scalarField& field) const
  {
    scalar qMin_=min(field);                                                                                        
    scalar qMax_=max(field);
    reduce(qMin_, minOp<scalar>());     
    reduce(qMin_, minOp<scalar>());
    Info<<this->patch().name()<<" "<<variableName<<":[ "<<qMin_<<"  "<<qMax_<<" ]"<<endl;       
  }
  void Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::minMax(word variableName,const scalarField& field) const
  {
    scalar qMin_=min(field);
    scalar qMax_=max(field);
    reduce(qMin_, minOp<scalar>());
    reduce(qMax_, maxOp<scalar>());                                                                 
    Info<<this->patch().name()<<" "<<variableName<<":[ "<<qMin_<<"  "<<qMax_<<" ]"<<endl;       
  }

void Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::write
(
 Ostream& os
 ) const
{
  mixedFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  namespace radiation
  {
    makePatchTypeField
    (
     fvPatchScalarField,
     wideBandDiffusiveRadiationMixedFvPatchScalarField
     );
  }
}


// ************************************************************************* //
