{
  vertT[0]=Twall_;
  // Calculate constant at edge of the surface Layer
  scalar zEdge_=0.2*PBLH;
  scalar phiM=(Foam::pos(MOL)*(1+ 5.0*zEdge_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*12*zEdge_/MOL,-1.0/3.0)));
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
  scalar RefC_=nutTemp*thetastar/(kappa_*zEdge_)*phiH;//Foam::pow(1-15*zEdge_/MOL,-1.0/2.0);
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
  scalar runCounter=0;
  if (solverType_=="steadyWind")
    lapseRate=1e-5;
  Info<<"Done Temperature:"<<endl;
  while(runCounter<10000)
    {
      forAll(zLoc,i)
	{
	  if(i<zLower)
	    vertT[i]=vertT[i];
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
	      vertT[i]=RefC_/alphaTemp*(zLoc[i]-zLoc[i-1])+vertT[i-1];
	    }
	  else if(i<zUpper)
	    {
	      scalar volume=0.5*(zLoc[i+1]-zLoc[i-1]);
	      scalar fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertT[i+1]-vertT[i])/(zLoc[i+1]-zLoc[i]);
	      scalar fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertT[i]-vertT[i-1])/(zLoc[i]-zLoc[i-1]);
	      vertT[i]=vertT[i]+0.01/(volume*Prt.value())*(fluxTop-fluxBottom);
	    }
	  else if(i==zUpper)
	    {
	      vertT[i]=TUpper;
	    }
	  else
	    vertT[i]=vertT[i-1];
	}
      runCounter=runCounter+1;
    }

  forAll(zLoc,i)
    {
      zProfile_[i]=zLoc[i];
      vertTProfile_[i]=vertT[i];
    }
  Info<<"Done Temperature:"<<endl;
}
