{
  //vertqV[0]=qvFace;
  // Compute qv at 
  //  x=Foam::pow(1- Foam::neg(MOL)*16*(0.01)/MOL,0.25);   
  //  psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*0.01/MOL);
  //  if(mag(MOL)>500)
  //    psiHL=0.0;
  //den=Foam::log(kappa_*ustar*0.01/ka+0.01/0.01)-psiHL;
  //  den=Foam::log(kappa_*ustar*0.01/ka+0.01/0.01)-psiHL;
  scalar phiH=(Foam::pos(MOL)*(1+ 5.0*zt_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*16*zt_/MOL,-1.0/2.0)));
  vertqV[0]=qvCell-qstar/kappa_*phiH;
  if(mAvail_>0.9)
	vertqV[0]=faceqsat;
  Info<<"Q*:"<<qstar<<"  "<<den<<"  "<<qvFace<<"  "<<vertqV[0]<<endl;
  // Calculate constant at edge of the surface Layer
  scalar zEdge_=0.2*PBLH;
  scalar phiM=(Foam::pos(MOL)*(1+ 5.0*zEdge_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*12*zEdge_/MOL,-1.0/3.0)));
  if(pos(MOL) && MOL<75)
    phiM=6;
  else if(mag(MOL)>500)
    phiM=1.0;
  scalar nutTemp=ustar*kappa_*zEdge_/phiM*Foam::sqr(1-zEdge_/PBLH);
  phiH=(Foam::pos(MOL)*(1+ 5.0*zEdge_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*15*zEdge_/MOL,-1.0/2.0)));
  if(pos(MOL) && MOL<75)
    phiH=6.0;
  else if(mag(MOL)>500)
    phiH=1.0;
  scalar RefC_=nutTemp*qstar/(kappa_*zEdge_)*phiH;//Foam::pow(1-15*zEdge_/MOL,-1.0/2.0);
  // Initialize Profile
  forAll(zLoc,i)
    {
      if(i>0)
	{
	  scalar zVal_=zLoc[i]+zt_;
	  phiH=(Foam::pos(MOL)*(1+ 5.0*zVal_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*16*zVal_/MOL,-1.0/2.0)));
	  vertqV[i]=vertqV[i-1]+qstar/kappa_*phiH;
	}
    }
  // Update T above surface layer 
  forAll(zLoc,i)
    {
      if(zLoc[i]>(0.2*PBLH))
	//if(i>0)
	{
	  scalar zTemp=min(zLoc[i],0.9*PBLH);
	  phiM=1;
	  if(pos(MOL) && MOL<75)
	    phiM=6;
	  else if(pos(MOL) && MOL<500)
	    phiM=1+5*zTemp/MOL;
	  else if(neg(MOL) && MOL>-500)
	    phiM=Foam::pow(1- 12*zTemp/MOL,-1.0/3.0);
	  nutTemp=ustar*kappa_*zTemp/phiM*Foam::sqr(1-zTemp/PBLH);
	  scalar alphaTemp=nutTemp/Prt.value();
	  if((zLoc[i]>(0.5*PBLH) && pos(MOL)) || vertT[i-1]>=TUpper)
	    vertqV[i]=vertqV[i-1];
	  else if ((zLoc[i]>(0.5*PBLH) && neg(MOL)) || vertT[i-1]<=TUpper)
	    vertqV[i]=vertqV[i-1];
	  else 
	    vertqV[i]=RefC_/alphaTemp*(zLoc[i]-zLoc[i-1])+vertqV[i-1];
	}
    }
  forAll(zLoc,i)
    {
      if(vertqV[i]<0)
	vertqV[i]=vertqV[i-1];
      zProfile_[i]=zLoc[i];
      vertqVProfile_[i]=vertqV[i];
      Info<<zProfile_[i]<<"    "<<vertTProfile_[i]<<" "<<vertqVProfile_[i]<<endl;
    }
}
