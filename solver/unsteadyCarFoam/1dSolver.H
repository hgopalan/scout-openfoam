// Compute u* from climate model data 
#include "ustar.H"
// Compute wall temperature and humidity 
#include "computeUstarThetastar.H"
// Initialize variables 
// nut 
scalar phiM=1.0;
forAll(zLoc,i)
{
  zProfile_[i]=zLoc[i];
  scalar z_ = zLoc[i]+Z0_;
  if(i==0)
    z_=Z0_;
  z_=min(z_,0.9*PBLH);
  scalar correction=Foam::sqr(1-Foam::min(z_,0.9*PBLH)/PBLH);
  if(pos(MOL) && MOL<500)
    phiM=1+5*z_/MOL;
  if(neg(MOL) && MOL>-500)
    phiM=Foam::pow(1- 12*z_/MOL,-1.0/3.0);
  if(pos(MOL) && MOL<2.5)
    vertNut[i]=ustar*kappa_*z_/6*correction;
  else
    vertNut[i]=ustar*kappa_*z_/phiM*correction;
  if(z_>(0.2*PBLH))
    vertNut[i]=max(vertNut[i],0.01);
}
// Wind 
forAll(zLoc,cellI)
{
  scalar z_ = zLoc[cellI]+Z0_;
  if(mag(MOL)>500)
    {
      vertU[cellI]=locFlowDir[0]*ustar/kappa_*(Foam::log(z_/Z0_)+z_/MBL-z_/PBLH*0.5*z_/MBL);
      vertV[cellI]=locFlowDir[1]*ustar/kappa_*(Foam::log(z_/Z0_)+z_/MBL-z_/PBLH*0.5*z_/MBL);
    }
  else if (pos(MOL) && MOL<2.5)
    {
      vertU[cellI]=locFlowDir[0]*ustar/kappa_*(6*Foam::log(z_/Z0_)-5*z_/PBLH+z_/MBL-z_/PBLH*0.5*z_/MBL);
      vertV[cellI]=locFlowDir[0]*ustar/kappa_*(6*Foam::log(z_/Z0_)-5*z_/PBLH+z_/MBL-z_/PBLH*0.5*z_/MBL);
    }
  else if(pos(MOL))
    {
      vertU[cellI]=locFlowDir[0]*ustar/kappa_*(Foam::log(z_/Z0_)+5*z_/MOL*(1-0.5*z_/PBLH)+z_/MBL-z_/PBLH*0.5*z_/MBL);
      vertV[cellI]=locFlowDir[1]*ustar/kappa_*(Foam::log(z_/Z0_)+5*z_/MOL*(1-0.5*z_/PBLH)+z_/MBL-z_/PBLH*0.5*z_/MBL);
    }
  else
    {
      scalar x=Foam::pow(1-12*z_/MOL,1.0/3.0);
      scalar psiMF=1.5*Foam::log((1+x+Foam::sqr(x))/3.0)+3.14/Foam::sqrt(3.0)-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0));
      scalar windSpeed=Foam::log(z_/Z0_)-psiMF+z_/PBLH*(1+(Foam::pow(1-12*z_/MOL,1.0/3.0)-1)/(8*z_/MOL));
      windSpeed=windSpeed+z_/MBL-z_/PBLH*0.5*z_/MBL;
      vertU[cellI]=locFlowDir[0]*ustar/kappa_*windSpeed;
      vertV[cellI]=locFlowDir[1]*ustar/kappa_*windSpeed;
    }
}
// Update i==0
vector UpParallel(vertU[1],vertV[1],0);
if(mag(MOL)>500)
  {
    scalar z_=2+Z0_;
    scalar deltaU=z_*ustar/kappa_*(1-z_/PBLH)*(1/z_+1/MBL+1/(PBLH-z_));
    UpParallel=(UpParallel*(1.0-deltaU/(mag(UpParallel)+SMALL)));
    vertU[0]=UpParallel[0];
    vertV[0]=UpParallel[1];
  }
 else if(pos(MOL) && MOL<2.5)
   {
     scalar z_=2+Z0_;
     scalar deltaU=z_*ustar/kappa_*(1-z_/PBLH)*(6/z_+1/MBL+1/(PBLH-z_));
     UpParallel=(UpParallel*(1.0-deltaU/(mag(UpParallel)+SMALL)));
     vertU[0]=UpParallel[0];
     vertV[0]=UpParallel[1];
   }
 else if(pos(MOL))
   {
     scalar z_=2+Z0_;
     scalar deltaU=z_*ustar/kappa_*(1-z_/PBLH)*((1+5*z_/MOL)/z_+1/MBL+1/(PBLH-z_));
     UpParallel=(UpParallel*(1.0-deltaU/(mag(UpParallel)+SMALL)));
     vertU[0]=UpParallel[0];
     vertV[0]=UpParallel[1];
   }
 else
   {
     scalar z_=2;
     scalar deltaU=ustar/kappa_*(z_+Z0_)/z_*Foam::pow(1-12*z_/MOL,-1.0/3.0);
     UpParallel=(UpParallel*(1.0-deltaU/(mag(UpParallel)+SMALL)));
     vertU[0]=UpParallel[0];
     vertV[0]=UpParallel[1];
   }
// Temperature 
vertT[0]=Twall_;
scalar zEdge_=0.2*PBLH;
phiM=(Foam::pos(MOL)*(1+ 5.0*zEdge_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*12*zEdge_/MOL,-1.0/3.0)));
if(pos(MOL) && MOL<2.5)
  phiM=6;
 else if(mag(MOL)>500)
   phiM=1.0;
scalar nutTemp=ustar*kappa_*zEdge_/phiM*Foam::sqr(1-zEdge_/PBLH);
scalar phiH=(Foam::pos(MOL)*(1+ 5.0*zEdge_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*15*zEdge_/MOL,-1.0/2.0)));
if(pos(MOL) && MOL<2.5)
  phiH=6.0;
 else if(mag(MOL)>500)
   phiH=1.0;
scalar RefT_=nutTemp*thetastar/(kappa_*zEdge_)*phiH;//Foam::pow(1-15*zEdge_/MOL,-1.0/2.0);
forAll(zLoc,i)
{
  scalar zTemp=zLoc[i]+zt_;	  
  if(zTemp> (0.2*PBLH) && lowerMarked==0)
    {
      zLower=i;
      lowerMarked=1;
    }
}
// Initialize Profile
forAll(zLoc,i)
{
  if(i>0)
    {
      scalar zTemp=min(zLoc[i],0.9*PBLH)+zt_;	  
      x=Foam::pow(1- Foam::neg(MOL)*16*(zTemp)/MOL,0.25);
      psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*zTemp/MOL);
      if(mag(MOL)>500)
	psiHL=0.0;
      vertT[i]=Twall_ + thetastar/kappa_*(Foam::log(zTemp/zt_)-psiHL);
      //	  if(i==1)
      if(pos(MOL) && MOL<2.5)
	vertT[i]=Twall_ + thetastar/kappa_*6*Foam::log(zTemp/zt_);
    }
}
// Vapor 
vertqV[0]=qwall_;
scalar RefqV_=nutTemp*qstar/(kappa_*zEdge_)*phiH;//Foam::pow(1-15*zEdge_/MOL,-1.0/2.0);
// Initialize Profile
forAll(zLoc,i)
{
  if(i>0)
    {
      scalar zTemp=min(zLoc[i],0.9*PBLH)+zt_;	  
      x=Foam::pow(1- Foam::neg(MOL)*16*(zTemp)/MOL,0.25);
      psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*zTemp/MOL);
      if(mag(MOL)>500)
	psiHL=0.0;
      vertqV[i]=qwall_ + qstar/kappa_*(Foam::log(zTemp/zt_)-psiHL);
      //	  if(i==1)
      if(pos(MOL) && MOL<2.5)
	vertqV[i]=qwall_ + qstar/kappa_*6*Foam::log(zTemp/zt_);
    }
}
forAll(zLoc,i)
{
  scalar zTemp=zLoc[i]+Z0_;
  if(zTemp> (1.1*PBLH) && upperMarked==0)
    {
      zUpper=i;
      upperMarked=1;
    }
  zTemp=zLoc[i]+zt_;	  
  if(zTemp> (0.2*PBLH) && lowerMarked==0)
    {
      zLower=i;
      lowerMarked=1;
    }        
}
scalar windCounter=0;
while(windCounter<50000)
  {
    forAll(zLoc,i)
      {
	scalar z_ = zLoc[i]+Z0_;
	/*if(MOL>0 && z_ >=(0.7*PBLH) && z_<(0.9*PBLH))
	  {
	    scalar invL=1/z_*(1+5*z_/MOL)+1/MBL+1/(PBLH-z_);
	    scalar l0_h=Foam::min(1/invL,50);
	    scalar MU=Foam::sqrt(Foam::sqr(vertU[i+1])+Foam::sqr(vertV[i+1]));
	    scalar ML=Foam::sqrt(Foam::sqr(vertU[i-1])+Foam::sqr(vertV[i-1]));
	    scalar dMdz=0.5*(MU-ML)/(zLoc[i+1]-zLoc[i-1]);
	    vertNut[i]=Foam::sqrt(Foam::sqr(l0_h)*Foam::mag(dMdz)*vertNutProfile_[i]);
	  }
	else if (MOL>0 && z_ >=(0.9*PBLH))
	vertNut[i]=vertNut[i-1];*/
	if(z_<(0.2*PBLH))
	  z_=z_;
	else if(i==zLower)
	  {
	    scalar zTemp=min(zLoc[i],0.9*PBLH)+zt_;
	    phiM=1;
	    if(pos(MOL) && MOL<2.5)
	      phiM=6;
	    else if(pos(MOL) && MOL<500)
	      phiM=1+5*zTemp/MOL;
	    else if(neg(MOL) && MOL>-500)
	      phiM=Foam::pow(1- 12*zTemp/MOL,-1.0/3.0);
	    nutTemp=ustar*kappa_*zTemp/phiM*Foam::sqr(1-zTemp/PBLH);
	    scalar alphaTemp=nutTemp/Prt.value();  
	    vertT[i]=RefT_/alphaTemp*(zLoc[i]-zLoc[i-1])+vertT[i-1];
	    vertqV[i]=RefqV_/alphaTemp*(zLoc[i]-zLoc[i-1])+vertqV[i-1];	      
	  }		
	else if(i>zLower && z_<(0.5*PBLH))
	  {
	    scalar volume=0.5*(zLoc[i+1]-zLoc[i-1]);
	    scalar fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertT[i+1]-vertT[i])/(zLoc[i+1]-zLoc[i]);
	    scalar fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertT[i]-vertT[i-1])/(zLoc[i]-zLoc[i-1]);		  
            vertT[i]=vertT[i]+0.01/(volume*Prt.value())*(fluxTop-fluxBottom);
	    volume=0.5*(zLoc[i+1]-zLoc[i-1]);
	    fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertqV[i+1]-vertqV[i])/(zLoc[i+1]-zLoc[i]);
	    fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertqV[i]-vertqV[i-1])/(zLoc[i]-zLoc[i-1]);
            vertqV[i]=vertqV[i]+0.01/(volume*Prt.value())*(fluxTop-fluxBottom);
	  }
	else if(z_ >=(0.5*PBLH) && z_<(1.1*PBLH))
	  {
	    scalar volume=0.5*(zLoc[i+1]-zLoc[i-1]);
	    scalar fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertU[i+1]-vertU[i])/(zLoc[i+1]-zLoc[i]);
	    scalar fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertU[i]-vertU[i-1])/(zLoc[i]-zLoc[i-1]);
	    vertU[i]=vertU[i]+0.01/volume*(fluxTop-fluxBottom);
	    fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertV[i+1]-vertV[i])/(zLoc[i+1]-zLoc[i]);
	    fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertV[i]-vertV[i-1])/(zLoc[i]-zLoc[i-1]);
	    vertV[i]=vertV[i]+0.01/volume*(fluxTop-fluxBottom);
	    volume=0.5*(zLoc[i+1]-zLoc[i-1]);
	    fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertT[i+1]-vertT[i])/(zLoc[i+1]-zLoc[i]);
	    fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertT[i]-vertT[i-1])/(zLoc[i]-zLoc[i-1]);		  
            vertT[i]=vertT[i]+0.01/(volume*Prt.value())*(fluxTop-fluxBottom);
	    volume=0.5*(zLoc[i+1]-zLoc[i-1]);
	    fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertqV[i+1]-vertqV[i])/(zLoc[i+1]-zLoc[i]);
	    fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertqV[i]-vertqV[i-1])/(zLoc[i]-zLoc[i-1]);
            vertqV[i]=vertqV[i]+0.01/(volume*Prt.value())*(fluxTop-fluxBottom);
	  }
	else if (i==zUpper)
	  {
	    vertU[i]=UUpper;
	    vertV[i]=VUpper;
	    vertT[i]=TUpper;
	    vertqV[i]=qVUpper;
	  }
	else
	  {
	    vertU[i]=vertU[i-1];
	    vertV[i]=vertV[i-1];
	    vertT[i]=vertT[i-1];
	    vertqV[i]=vertqV[i-1];
	  }
      }
    windCounter=windCounter+1;
  }
forAll(zLoc,i)
{
  vertNutProfile_[i]=vertNut[i];
  vertUProfile_[i]=vertU[i];
  vertVProfile_[i]=vertV[i];
  vertTProfile_[i]=vertT[i];
  vertqVProfile_[i]=vertqV[i];
  //  if(zLoc[i]<1.2*PBLH)
    //    Info<<i<<"  "<<zProfile_[i]<<"    "<<vertUProfile_[i]<<"  "<<vertVProfile_[i]<<"   "
    //	<<vertNutProfile_[i]<<"  "<<vertTProfile_[i]<<" "<<vertqVProfile_[i]<<endl;	  
}     
