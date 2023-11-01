{
  // Compute Wind 
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
  forAll(zLoc,i)
    {
      vertUProfile_[i]=vertU[i];
      vertVProfile_[i]=vertV[i];
    }

  forAll(zLoc,i)
    {
      scalar zTemp=zLoc[i]+Z0_;
      if(zTemp> (1.1*PBLH) && upperMarked==0)
        {
          zUpper=i;
          upperMarked=1;
        }
    }
  Info<<"Upper:"<<zUpper<<endl;
  if (solverType_=="steadyWind")
    lapseRate=0;
  scalar windCounter=0;
  while(windCounter<10000)
    {
      forAll(zLoc,i)
	{
	  scalar z_ = zLoc[i]+Z0_;
	  if(z_<(0.5*PBLH))
	    z_=z_;      
	  else if(z_ >=(0.5*PBLH) and z_<(1.1*PBLH))
	    {
	      scalar volume=0.5*(zLoc[i+1]-zLoc[i-1]);
	      scalar fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertU[i+1]-vertU[i])/(zLoc[i+1]-zLoc[i]);
	      scalar fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertU[i]-vertU[i-1])/(zLoc[i]-zLoc[i-1]);
	      vertU[i]=vertU[i]+0.1/volume*(fluxTop-fluxBottom);
	      fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertV[i+1]-vertV[i])/(zLoc[i+1]-zLoc[i]);
	      fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertV[i]-vertV[i-1])/(zLoc[i]-zLoc[i-1]);
	      vertV[i]=vertV[i]+0.1/volume*(fluxTop-fluxBottom);
	    }
	  else if (i==zUpper)
            {
              vertU[i]=UUpper;
              vertV[i]=VUpper;
            }
	  else
	    {
	      vertU[i]=vertU[i-1];
	      vertV[i]=vertV[i-1];
	    }
	}
      windCounter=windCounter+1;
    }
  forAll(zLoc,i)
    {
      vertUProfile_[i]=vertU[i];
      vertVProfile_[i]=vertV[i];
    }     
}

