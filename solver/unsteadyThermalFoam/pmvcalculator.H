{
  // PMV Calculator - Based on CBE Thermal Tool 
  //const volScalarField swr_=U.db().lookupObject<volScalarField>("QrSolar"); 
  //const volScalarField lwr_=U.db().lookupObject<volScalarField>("qr"); 
  //scalar sigmaSB= 5.670367e-8;
  //surfaceScalarField swrsurface_=fvc::interpolate(swr_);
  //surfaceScalarField lwrsurface_=fvc::interpolate(lwr_);
  //surfaceScalarField Tmrtsurface_=pow(max(swrsurface_+lwrsurface_,0*lwrsurface_)/sigmaSB,0.25)-273.15;
  //Tmrt=fvc::average(Tmrtsurface_);
  //forAll(Tmrt,cellI)
  //{
  //if(Tmrt[cellI]<T[cellI])
  //Tmrt[cellI]=T[cellI];
  //}
  //  Tmrt=thetaV;
  const volScalarField swr_=U.db().lookupObject<volScalarField>("qmrtsw");
  const volScalarField lwr_=U.db().lookupObject<volScalarField>("qmrtlw");
  scalar nightCorr=1;
  if(MOL>0)
    nightCorr=0;
  scalar sigmaSB= 5.670367e-8; 
  volScalarField radiation_=max(nightCorr*swr_+lwr_,0*lwr_);
  forAll(radiation_,cellI)
    {
      if(radiation_[cellI]>1200)
	radiation_[cellI]=1200;
      if(radiation_[cellI]<-400)
	radiation_[cellI]=-400;
    }
  radiation_.correctBoundaryConditions();
  Tmrt=pow(radiation_/sigmaSB,0.25);//-273.15;
  Tmrt=max(Tmrt,thetaV);
  Tmrt.correctBoundaryConditions();
  scalar wme = 0;
  // CBE Inputs 
  scalar clo = 0.61;
  scalar met = 1.2;
  //  scalar rh = 80;
  //
  forAll(Tmrt,cellI)
    {
      scalar ta=T[cellI]-273.15;
      scalar tr=Tmrt[cellI]-273.15;
      scalar vel=mag(U[cellI]);
      scalar pa = relHum[cellI] * 10 * Foam::exp(16.6536 - 4030.183 / (ta + 235));
      scalar icl = 0.155 * clo;
      scalar m = met * 58.15 ;
      scalar w = wme * 58.15;  
      scalar mw = m - w;
      scalar fcl= 1.05 + (0.645 * icl);
      if (icl <= 0.078)
        fcl = 1 + (1.29 * icl);
      scalar hcf = 12.1 * Foam::sqrt(vel);
      scalar taa = ta + 273;
      scalar tra=tr + 273;
      scalar tcla = taa + (35.5 - ta) / (3.5 * icl + 0.1);
      scalar p1 = icl * fcl;
      scalar p2 = p1 * 3.96;
      scalar p3 = p1 * 100;
      scalar p4 = p1 * taa;
      scalar p5 = (308.7 - 0.028 * mw) + (p2 * Foam::pow(tra / 100, 4));
      scalar xn = tcla / 100.0;
      scalar xf = tcla / 50.0;
      scalar eps = 0.00015;
      scalar n = 0;
      scalar hc=0;
      scalar hcn=0;
      while (mag(xn - xf) > eps)
	{
	  xf = (xf + xn) / 2;
	  hcn = 2.38 * Foam::pow(mag(100.0 * xf - taa), 0.25);
	  if (hcf > hcn)
            hc = hcf;
	  else
            hc = hcn;
	  xn = (p5 + p4 * hc - p2 * Foam::pow(xf, 4)) / (100 + p3 * hc);
	  n=n+1;
	  if (n > 150)
	    {
	      Info<<"Iterations exceeded:"<<endl;
	      break;
	    }
	}

      scalar tcl = 100 * xn - 273;
      scalar hl1 = 3.05 * 0.001 * (5733 - (6.99 * mw) - pa);
      scalar hl2 = 0.0;
      if (mw > 58.15)
        hl2 = 0.42 * (mw - 58.15);
      scalar hl3 = 1.7 * 0.00001 * m * (5867 - pa);
      scalar hl4 = 0.0014 * m * (34 - ta);
      scalar hl5 = 3.96 * fcl * (Foam::pow(xn, 4) - Foam::pow(tra / 100, 4));
      scalar hl6 = fcl * hc * (tcl - ta);
      scalar ts = 0.303 * Foam::exp(-0.036 * m) + 0.028;
      PMV[cellI]=min(ts * (mw - hl1 - hl2 - hl3 - hl4 - hl5 - hl6),3);
      //ppd = 100.0 - 95.0 * math.exp(-0.03353 * pow(pmv, 4.0)
      //- 0.2179 * pow(pmv, 2.0))
    }
  bound(PMV,-3);
  PMV.correctBoundaryConditions();
  Info<<"MRT Bounds:"<<gMin(Tmrt)<<"    "<<gMax(Tmrt)<<endl;   
  Info<<"PMV Bounds:"<<gMin(PMV)<<"    "<<gMax(PMV)<<endl; 
}
