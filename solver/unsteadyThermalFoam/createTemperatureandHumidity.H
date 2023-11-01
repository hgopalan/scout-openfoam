if(mesh.time().value()==0 || (word(U.mesh().ddtScheme("default")) == "steadyState") || generateThermal==1)
  {
    forAll(T,cellI)
      {
	vector cc=mesh.C()[cellI];
	scalar ccz=cc[2];
	T[cellI]=interpolateXY(ccz,zLoc,vertT);
      }
    //    T.write();
    // Humidity
    /*scalar cellesat=610.94*Foam::exp(17.625*(TRef-273.15)/(TRef-30.11));
      scalar cellqsat=(0.622 *cellesat/(101325-0.378*cellesat));
      scalar faceesat=610.94*Foam::exp(17.625*(Twall_-273.15)/(Twall_-30.11));
      scalar faceqsat=(0.622 *faceesat/(101325-0.378*faceesat));
      scalar qvCell=RH*cellqsat;                                                                                                  
      scalar qvFace=RH*faceqsat;
      x=Foam::pow(1- Foam::neg(MOL)*16*(TZRef_)/MOL,0.25);                                                                        
      psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*TZRef_/MOL);                                                     
      scalar qstar=(qvCell-qvFace)*kappa_/(Foam::log(TZRef_+surfacez0/surfacez0)-psiHL);       */
    /*forAll(relHum,cellI)
      {
      scalar z_=y_[cellI];
      if(z_<=(0.2*PBLH))
      {
      x=Foam::pow(1- Foam::neg(MOL)*16*(z_)/MOL,0.25);
      psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*z_/MOL);
      qv[cellI]=qvFace+qstar/kappa_*(Foam::log((z_+ Z0_)/Z0_)-psiHL);
      }
      else
      {
      z_=0.2*PBLH;
      x=Foam::pow(1- Foam::neg(MOL)*16*(z_)/MOL,0.25);
      psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*z_/MOL);
      qv[cellI]=qvFace+qstar/kappa_*(Foam::log((z_+ Z0_)/Z0_)-psiHL);
      }
      scalar esat=610.94*Foam::exp(17.625*(T[cellI]-273.15)/(T[cellI]-30.11));                                                   
      scalar qsat=(0.622*esat/(101325-0.378*esat));
      relHum[cellI]=qv[cellI]/qsat;     
      }*/

    forAll(qv,cellI)
      {
	vector cc=mesh.C()[cellI];
	scalar ccz=cc[2];
	qv[cellI]=interpolateXY(ccz,zLoc,vertqV);
	scalar esat=610.94*Foam::exp(17.625*(T[cellI]-273.15)/(T[cellI]-30.11));
	scalar qsat=(0.622*esat/(101325-0.378*esat));  
	relHum[cellI]=qv[cellI]/qsat;     
      }
    qv.correctBoundaryConditions();
    //    qv.write();
    relHum.correctBoundaryConditions();
    //    relHum.write();
  }
T.correctBoundaryConditions();//We do not want to screw up BC !!
thetaV=T*(1+(461.51/287-1)*qv);
//scalar cellesat=610.94*Foam::exp(17.625*(TRef-273.15)/(TRef-30.11));
//scalar cellqsat=(0.622 *cellesat/(101325-0.378*cellesat));
scalar qvRef=RH;
TRef=TRef*(1+(461.51/287-1)*qvRef);
volScalarField volTRef=T*(1+(461.51/287-1)*qvRef);
rhok= 1 -1.0/volTRef*(T*(1+(461.51/287-1)*qv)-volTRef);
TRef=interpolateXY(mesh.time().value(),timevalues,TRefvalues);
//rhok= 1 -1.0/TRef*(thetaV-TRef);
// Temperature
Info<<"Potential Temperature Bounds:"<<gMin(T)<<"    "<<gMax(T)<<endl;
Info<<"Virtual Potential Temperature Bounds:"<<gMin(thetaV)<<"    "<<gMax(thetaV)<<endl;
// Specific Humidity
Info<<" Specific Humidity Bounds:"<<gMin(qv)<<"  "<<gMax(qv)<<endl;
Info<<" Relative  Humidity Bounds:"<<gMin(relHum)<<"  "<<gMax(relHum)<<endl;
generateThermal=0;
generateWind=0;



