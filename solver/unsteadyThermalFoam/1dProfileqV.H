{
  vertqV[0]=qwall_;
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
  scalar RefC_=nutTemp*qstar/(kappa_*zEdge_)*phiH;//Foam::pow(1-15*zEdge_/MOL,-1.0/2.0);
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
  scalar runCounter=0;
  if (solverType_=="steadyWind")
    lapseRate=2e-8;
  while(runCounter<10000)
    {
      forAll(zLoc,i)
	{
	  if(i<zLower)
	    vertqV[i]=vertqV[i];
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
	      vertqV[i]=RefC_/alphaTemp*(zLoc[i]-zLoc[i-1])+vertqV[i-1];
	    }
	  else if (i<zUpper)
	    {
	      scalar volume=0.5*(zLoc[i+1]-zLoc[i-1]);
	      scalar fluxTop=0.5*(vertNut[i+1]+vertNut[i])*(vertqV[i+1]-vertqV[i])/(zLoc[i+1]-zLoc[i]);
	      scalar fluxBottom=0.5*(vertNut[i]+vertNut[i-1])*(vertqV[i]-vertqV[i-1])/(zLoc[i]-zLoc[i-1]);
	      vertqV[i]=vertqV[i]+0.01/(volume*Prt.value())*(fluxTop-fluxBottom);
	    }
	  else if(i==zUpper)
	    vertqV[i]=qVUpper;
	  else
	    vertqV[i]=vertqV[i-1];
	}
      runCounter=runCounter+1;
    }

  forAll(zLoc,i)
    {
      zProfile_[i]=zLoc[i];
      vertqVProfile_[i]=vertqV[i];
      if(zLoc[i]<1.5*PBLH)
	Info<<i<<"  "<<zProfile_[i]<<"    "<<vertUProfile_[i]<<"  "<<vertVProfile_[i]<<"   "
	    <<vertNutProfile_[i]<<"  "<<vertTProfile_[i]<<" "<<vertqVProfile_[i]<<endl;
    }
}
