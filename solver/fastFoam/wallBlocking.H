// Cavity Correction
Info<<"Correcting building zone"<<endl;
// Wake Correction
forAll(buildingZone,cellI)
{
  if(buildingZone[cellI]>=0)
    {
      scalar id=buildingZone[cellI];
      scalar xmin=buildingminX[id];
      scalar xmax=buildingmaxX[id];
      scalar ymin=buildingminY[id];
      scalar ymax=buildingmaxY[id];
      scalar zmin=buildingminZ[id];
      scalar zmax=buildingmaxZ[id];
      // Kaplan and Dinar 
      scalar L=xmax-xmin;
      scalar W=ymax-ymin;
      scalar H=zmax-zmin;
      vector cc=mesh.C()[cellI];
      scalar Lf=H*2*(W/H)/(1+0.8*W/H);
      scalar Lr=H*1.8*W/H/(Foam::pow(L/H,0.3)*(1+0.24*W/H));
      // Check for displacement zone 
      if(cc[0]<xmin && cc[0]>(xmin-Lf))
	{
	  scalar X=xmin-cc[0];
	  scalar Y=0.5*(ymin+ymax)-cc[1];
	  scalar Z=0.5*(zmin+zmax)-cc[2];
	  scalar refvalue=Foam::sqr(Y/W)+Foam::sqr(X/Lf)*1/(1-Foam::sqr(Z/(0.6*H)));
	  if(refvalue>=0 && refvalue<=1)
	    {
	      U[cellI]=0*locFlowDir;
	      blankcell[cellI]=refvalue;
	    }
	}
      if(cc[0]>xmax && cc[0]<(xmax+Lr))
	{
	  scalar X=cc[0]-xmax;
	  scalar Y=0.5*(ymin+ymax)-cc[1];
	  scalar Z=0.5*(zmin+zmax)-cc[2];
	  scalar refvalue=Foam::sqr(Y/W)+Foam::sqr(X/Lr)*1/(1-Foam::sqr(Z/H));
	  scalar dN=Lr*Foam::sqrt((1-Foam::sqr(Z/H))*(1-Foam::sqr(Y/W)))-0.5*L;
	  scalar x=Foam::pow(1-Foam::neg(MOL)*12*ZREF_/MOL,1.0/3.0);
	  scalar psiM=(Foam::pos(MOL)*(-5.0* ZREF_/MOL)+Foam::neg(MOL)*(1.5*Foam::log(1.0/3.0*(1+x+Foam::sqr(x)))-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0))+3.14/Foam::sqrt(3.0)));
	  scalar ustar= M*kappa_/(Foam::log((ZREF_+Z0_)/Z0_)-psiM+SMALL);
	  // Roughness correction
	  ustar=ustar*Foam::log(260/Z0_)/Foam::log(260/surfacez0);
	  // Calculate wind speed 
	  scalar z_ = H;
	  z_=min(z_,0.9*PBLH)+surfacez0;
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
	  vector urefwind = upperstablecorrection*locFlowDir*ustar/kappa_*(Foam::log((z_+surfacez0)/surfacez0)-psiM+extraterm+SMALL);
	  if(refvalue>=0 && refvalue<=1)
	    {
	      U[cellI]=-urefwind*Foam::sqr(1-X/dN);
	      blankcell[cellI]=refvalue;
	    }
	}
      if(cc[0]>=(xmax+Lr) && cc[0]<(xmax+4*Lr))
	{
	  scalar X=cc[0]-xmax;
	  scalar Y=0.5*(ymin+ymax)-cc[1];
	  scalar Z=0.5*(zmin+zmax)-cc[2];
	  scalar refvalue=Foam::sqr(Y/W)+Foam::sqr(X/(3*Lr))*1/(1-Foam::sqr(Z/H));
	  scalar dN=Lr*Foam::sqrt((1-Foam::sqr(Z/H))*(1-Foam::sqr(Y/W)))-0.5*L;
	  if(refvalue>=0 && refvalue<=1)
	    {
	      blankcell[cellI]=refvalue;
	      U[cellI]*=(1-Foam::pow(dN/X,1.5));
	    }
	}
      scalar R=Foam::pow(Foam::min(W,H),0.67)*Foam::pow(Foam::min(W,H),0.33);
      scalar Ht=0.22*R;
      scalar Lt=0.9*R;
      if((cc[0]>xmin && cc[0]<xmax) && (cc[1]>ymin && cc[1]<ymax) && cc[2]>(0.9*zmax))
	{ 
	  scalar X=0.5*(xmin+xmax)-cc[0];
	  scalar Y=0.5*(ymin+ymax)-cc[1];
	  scalar refvalue=Foam::sqr(X/L)+Foam::sqr(Y/W);
	  if ((cc[2]-zmax)<0.5*Ht && (cc[0]-xmin)<Lt && refvalue>=0 && refvalue<=1)
	    {
	      U[cellI][0]=-U[cellI][0]*(1+0.05*W/H);
	    }
	  else if((cc[2]-zmax)>0.5*Ht && (cc[0]-xmin)<Lt && refvalue>=0 && refvalue<=1)
	    {
	      U[cellI][0]=U[cellI][0]*(1+0.05*W/H);		 
	    }
	}

    }
} 
U.correctBoundaryConditions();
blankcell.write();

