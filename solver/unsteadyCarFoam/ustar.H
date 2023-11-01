scalar denominator=Foam::log(ZREF_/Z0_)+ZREF_/MBL-ZREF_/PBLH*(0.5*ZREF_/MBL);
scalar ustar=M*kappa_/denominator;
if(pos(MOL) && MOL>=2.5 && MOL<500)
  {
    scalar denominator=Foam::log(ZREF_/Z0_)+5*ZREF_/MOL*(1-0.5*ZREF_/PBLH)+ZREF_/MBL-ZREF_/PBLH*0.5*ZREF_/MBL;
    ustar= M*kappa_/denominator;
  }
 else if(pos(MOL) && MOL<2.5)
   {
     scalar denominator=6*Foam::log(ZREF_/Z0_)-5*ZREF_/PBLH+ZREF_/MBL-ZREF_/PBLH*0.5*ZREF_/MBL;
     ustar=M*kappa_/denominator;
     //     Info<<"Stable Ustar:"<<ustar<<endl;
   }
 else
   {
     scalar x=Foam::pow(1-12*ZREF_/MOL,1.0/3.0);
     scalar psiM=1.5*Foam::log((1+x+Foam::sqr(x))/3.0)+3.14/Foam::sqrt(3.0)-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0));
     scalar denominator=Foam::log(ZREF_/Z0_)-psiM+ZREF_/PBLH*(1+(Foam::pow(1-12*ZREF_/MOL,1.0/3.0)-1)/(8*ZREF_/MOL));
     denominator=denominator+ZREF_/MBL-ZREF_/PBLH*0.5*ZREF_/MBL;
     ustar=M*kappa_/denominator;
   }
//ustar=max(ustar,0.1);
