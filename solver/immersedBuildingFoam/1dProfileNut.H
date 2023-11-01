scalar phiM=1.0;
{
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
      if(i>1)
	vertNut[i]=max(vertNut[i],0.1);
      vertNutProfile_[i]=vertNut[i];
    }
}
