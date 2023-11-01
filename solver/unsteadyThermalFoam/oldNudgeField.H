{
  scalar kappa_ = 0.41;
  scalar locMOL=interpolateXY(U.time().value(),timevalues,movalues);
  scalar locU1=interpolateXY(U.time().value(),timevalues,uvalues);
  scalar locV1=interpolateXY(U.time().value(),timevalues,vvalues);
  scalar locT=interpolateXY(U.time().value(),timevalues,TRefvalues);
  scalar locM=Foam::sqrt(sqr(locU1)+sqr(locV1));
  scalar xZ_=1.0;
  if(locMOL<0)
    xZ_ = Foam::pow(1- 16*ZREF_/locMOL,0.25);
  scalar intermediateZ_ = (locMOL < 0) ? (Foam::log(0.5*(1+sqr(xZ_))*sqr(0.5*(1+xZ_)))-2*Foam::atan(xZ_)+0.5*3.14) : (-5.0* ZREF_/locMOL) ; 
  vector locFlowDir(locU1/locM,locV1/locM,0);
  scalar locUstar=locM*kappa_/(Foam::log((ZREF_+Z0_)/Z0_)-intermediateZ_+SMALL); 
  if(locMOL<0)
    xZ_ = Foam::pow(1- 16*TZRef_/locMOL,0.25);
  scalar intermediateZT_ = (locMOL < 0) ? (2*Foam::log(0.5*(1+Foam::sqr(xZ_)))) : (-5.0* TZRef_/locMOL) ;
  forAll(verticalzprofile_,i)
    {
      scalar z_ = verticalzprofile_[i];
      scalar x_=0.0;
      if(locMOL<0)
	x_ = Foam::pow(1- 16*z_/locMOL,0.25);
      scalar intermediateU_ = (locMOL < 0) ? (Foam::log(0.5*(1+sqr(x_))*sqr(0.5*(1+x_)))-2*Foam::atan(x_)+0.5*3.14) : (-5.0* z_/locMOL) ;
      verticaluprofile_[i]=locFlowDir[0]*locUstar/kappa_*(Foam::log((z_+Z0_)/Z0_)-intermediateU_+SMALL);
      verticalvprofile_[i]=locFlowDir[1]*locUstar/kappa_*(Foam::log((z_+Z0_)/Z0_)-intermediateU_+SMALL);
      scalar intermediateT_ = (locMOL < 0) ? (2*Foam::log(0.5*(1+Foam::sqr(x_)))) : (-5.0* z_/locMOL) ;
      scalar thetastar=locT*Foam::sqr(locUstar)/(kappa_*9.81*locMOL);
      scalar Twall_= locT - thetastar/kappa_*(Foam::log((TZRef_+Z0_)/Z0_)-intermediateZT_);
      verticalTprofile_[i]=Twall_+thetastar/kappa_*(Foam::log((z_+ Z0_)/Z0_)-intermediateT_);
    }
  forAll(verticalzprofile_,i)
    {
      verticalzvalues[i] = verticalzprofile_[i];
      verticalURefvalues[i]=verticaluprofile_[i];
      verticalVRefvalues[i]=verticalvprofile_[i];
      verticalTRefvalues[i]=verticalTprofile_[i];
    }
}
#include "computeNewAverageFields.H"
forAll(U,cellI)
{
  scalar zc=mesh.C()[cellI][2];
  //U_avg[cellI][0]=U[cellI][0];
  //U_avg[cellI][1]=U[cellI][1];
  //T_avg[cellI]=T[cellI];
  U_avg[cellI][0]=interpolateXY(zc,verticalzvalues,uavg);
  U_avg[cellI][1]=interpolateXY(zc,verticalzvalues,vavg);
  Unudge[cellI][0]=interpolateXY(zc,verticalzvalues,verticalURefvalues);
  Unudge[cellI][1]=interpolateXY(zc,verticalzvalues,verticalVRefvalues);
}

