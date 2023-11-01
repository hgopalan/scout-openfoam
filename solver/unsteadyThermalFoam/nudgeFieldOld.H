{
scalar x=Foam::pow(1-Foam::neg(MOL)*12*ZREF_/MOL,1.0/3.0);
scalar psiM=(Foam::pos(MOL)*(-5.0* ZREF_/MOL)+Foam::neg(MOL)*(1.5*Foam::log(1.0/3.0*(1+x+Foam::sqr(x)))-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0))+3.14/Foam::sqrt(3.0)));
scalar ustar= M*kappa_/(Foam::log((ZREF_+Z0_)/Z0_)-psiM+SMALL);
scalar thetastar=TRef*Foam::sqr(ustar)/(kappa_*9.81*MOL);
x=Foam::pow(1- Foam::neg(MOL)*16*(TZRef_)/MOL,0.25);
scalar psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*TZRef_/MOL);
scalar Twall_= TRef - thetastar/kappa_*(Foam::log((TZRef_+Z0_)/Z0_)-psiHL);
scalar tSky=0.0555*Foam::pow(TRef,1.5);
scalar zSky=2500.0;
scalar zSurfaceTop=0.2*PBLH;
// Compute Wind 
forAll(verticalzprofile_,cellI)
{
  // Calculate wind speed 
  scalar z_ = verticalzprofile_[cellI];
  z_=min(z_,0.9*PBLH)+Z0_;
  x=Foam::pow(1-Foam::neg(MOL)*12*z_/MOL,1.0/3.0);
  psiM=(Foam::pos(MOL)*(-5.0* z_/MOL)+Foam::neg(MOL)*(1.5*Foam::log(1.0/3.0*(1+x+Foam::sqr(x)))-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0))+3.14/Foam::sqrt(3.0)));
  scalar extraterm=0;
  scalar upperstablecorrection=1;
  if(mag(MOL)>500)
    extraterm=z_/MBL-z_/PBLH*0.5*z_/MBL;    
  else if(pos(MOL))
    {
      if(y_[cellI]>PBLH)
	upperstablecorrection=Foam::min(y_[cellI]/(0.9*PBLH),1.2);
      extraterm=z_/MBL-z_/PBLH*0.5*z_/MBL-5*z_/MOL*0.5*z_/PBLH;    
    }
  else
    {
      extraterm=z_/MBL-z_/PBLH*0.5*z_/MBL+z_/PBLH*(1+(Foam::pow(1-12*z_/MOL,2.0/3.0)-1)/(8*z_/MOL));
    }  
  vector verticalU = upperstablecorrection*locFlowDir*ustar/kappa_*(Foam::log((z_+Z0_)/Z0_)-psiM+extraterm+SMALL);
  verticaluprofile_[cellI]=verticalU[0];
  verticalvprofile_[cellI]=verticalU[1];
}
  forAll(verticalzprofile_,i)
    {
      verticalzvalues[i] = verticalzprofile_[i];
      verticalURefvalues[i]=verticaluprofile_[i];
      verticalVRefvalues[i]=verticalvprofile_[i];
      //verticalTRefvalues[i]=verticalTprofile_[i];
    }
}
#include "computeNewAverageFields.H"
forAll(U,cellI)
{
  scalar zc=mesh.C()[cellI][2];
  U_avg[cellI][0]=interpolateXY(zc,verticalzvalues,uavg);
  U_avg[cellI][1]=interpolateXY(zc,verticalzvalues,vavg);
  //  Unudge[cellI][0]=interpolateXY(zc,verticalzvalues,verticalURefvalues);
  //  Unudge[cellI][1]=interpolateXY(zc,verticalzvalues,verticalVRefvalues);
}

