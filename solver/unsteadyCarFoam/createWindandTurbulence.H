Info<<"Creating Initial Conditions"<<endl;
// Computer ustar
//#include "ustar.H"
//#include "1dProfileNut.H"
//#include "1dProfileWind.H"
/*scalar denominator=Foam::log(ZREF_/Z0_)+ZREF_/MBL-ZREF_/PBLH*(0.5*ZREF_/MBL);
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
     Info<<"Stable Ustar:"<<ustar<<endl;
   }
 else 
   {
     scalar x=Foam::pow(1-12*ZREF_/MOL,1.0/3.0);
     scalar psiM=1.5*Foam::log((1+x+Foam::sqr(x))/3.0)+3.14/Foam::sqrt(3.0)-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0));
     scalar denominator=Foam::log(ZREF_/Z0_)-psiM+ZREF_/PBLH*(1+(Foam::pow(1-12*ZREF_/MOL,1.0/3.0)-1)/(8*ZREF_/MOL));
     denominator=denominator+ZREF_/MBL-ZREF_/PBLH*0.5*ZREF_/MBL;
     ustar=M*kappa_/denominator;
   } */
volVectorField Unudge=U;
volScalarField nutnudge=nut;
if(mesh.time().value()==0 || (word(U.mesh().ddtScheme("default")) == "steadyState") || generatewind==1)
  {
    // Compute Wind 
    /*    forAll(mesh.C(),cellI)
      {
	vector cc=mesh.C()[cellI];
	scalar z_ = cc[2]+Z0_;
	//  scalar z_=y_[cellI];
	z_=min(z_,0.9*PBLH);
	if(mag(MOL)>500)
	  U[cellI]=locFlowDir*ustar/kappa_*(Foam::log(z_/Z0_)+z_/MBL-z_/PBLH*0.5*z_/MBL);
	else if (pos(MOL) && MOL<2.5)
	  {
	    U[cellI]=locFlowDir*ustar/kappa_*(6*Foam::log(z_/Z0_)-5*z_/PBLH+z_/MBL-z_/PBLH*0.5*z_/MBL);
	  }
	else if(pos(MOL))
	  {
	    U[cellI]=locFlowDir*ustar/kappa_*(Foam::log(z_/Z0_)+5*z_/MOL*(1-0.5*z_/PBLH)+z_/MBL-z_/PBLH*0.5*z_/MBL);
	  }
	else
	  {
	    scalar x=Foam::pow(1-12*z_/MOL,1.0/3.0);
	    scalar psiMF=1.5*Foam::log((1+x+Foam::sqr(x))/3.0)+3.14/Foam::sqrt(3.0)-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0));
	    scalar windSpeed=Foam::log(z_/Z0_)-psiMF+z_/PBLH*(1+(Foam::pow(1-12*z_/MOL,1.0/3.0)-1)/(8*z_/MOL));
	    windSpeed=windSpeed+z_/MBL-z_/PBLH*0.5*z_/MBL;
	    U[cellI]=locFlowDir*ustar/kappa_*windSpeed;
	  }
	  } */
    forAll(U,cellI)
      {
        vector cc=mesh.C()[cellI];
        scalar ccz=cc[2];
        U[cellI][0]=interpolateXY(ccz,zLoc,vertU);
        U[cellI][1]=interpolateXY(ccz,zLoc,vertV);
	U[cellI][2]=0.0;
      }
    U.correctBoundaryConditions();
    //    U.write();
    // Compute Turbulence
    forAll(nut,cellI)
      {
        vector cc=mesh.C()[cellI];
        scalar ccz=cc[2];
        nut[cellI]=interpolateXY(ccz,zLoc,vertNut);
	alphaT[cellI]=nut[cellI]/0.7;
      }
    nut.correctBoundaryConditions();
    //    nut.write();
    alphaT.correctBoundaryConditions();
    //    alphaT.write();
    Unudge=U;
    nutnudge=nut;
    phi=fvc::flux(U);
    p_rgh=0*p_rgh;
  }
 else
   {
     // Compute Wind 
     /*forAll(mesh.C(),cellI)
       {
	 vector cc=mesh.C()[cellI];
	 scalar z_ = cc[2]+Z0_;
	 //  scalar z_=y_[cellI];
	 z_=min(z_,0.9*PBLH);
	 if(mag(MOL)>500)
	   Unudge[cellI]=locFlowDir*ustar/kappa_*(Foam::log(z_/Z0_)+z_/MBL-z_/PBLH*0.5*z_/MBL);
	 else if (pos(MOL) && MOL<2.5)
	   {
	     Unudge[cellI]=locFlowDir*ustar/kappa_*(6*Foam::log(z_/Z0_)-5*z_/PBLH+z_/MBL-z_/PBLH*0.5*z_/MBL);
	   }
	 else if(pos(MOL))
	   {
	     Unudge[cellI]=locFlowDir*ustar/kappa_*(Foam::log(z_/Z0_)+5*z_/MOL*(1-0.5*z_/PBLH)+z_/MBL-z_/PBLH*0.5*z_/MBL);
	   }
	 else
	   {
	     scalar x=Foam::pow(1-12*z_/MOL,1.0/3.0);
	     scalar psiMF=1.5*Foam::log((1+x+Foam::sqr(x))/3.0)+3.14/Foam::sqrt(3.0)-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0));
	     scalar windSpeed=Foam::log(z_/Z0_)-psiMF+z_/PBLH*(1+(Foam::pow(1-12*z_/MOL,1.0/3.0)-1)/(8*z_/MOL));
	     windSpeed=windSpeed+z_/MBL-z_/PBLH*0.5*z_/MBL;
	     Unudge[cellI]=locFlowDir*ustar/kappa_*windSpeed;
	   }
	   } */
     forAll(U,cellI)
       {
	 vector cc=mesh.C()[cellI];
	 scalar ccz=cc[2];
	 Unudge[cellI][0]=interpolateXY(ccz,zLoc,vertU);
	 Unudge[cellI][1]=interpolateXY(ccz,zLoc,vertV);
	 Unudge[cellI][2]=0.0;
       }
     // Compute Turbulence
     forAll(nut,cellI)
       {
	 vector cc=mesh.C()[cellI];
	 scalar ccz=cc[2];
	 nutnudge[cellI]=interpolateXY(ccz,zLoc,vertNut);
       }
   }	
generatewind=0;
Info<<"Wind Bounds:"<<gMin(U)<<"    "<<gMax(U)<<endl;
Info<<"Turbulence Bounds:"<<gMin(nut)<<"    "<<gMax(nut)<<endl;
