Info<<"TRef:"<<TRef<<endl;
scalar kappa_ = 0.41;
scalar locU1=interpolateXY(U.time().value(),timevalues,uvalues);
scalar locV1=interpolateXY(U.time().value(),timevalues,vvalues);
scalar M=Foam::sqrt(sqr(locU1)+sqr(locV1));
vector locFlowDir(locU1/M,locV1/M,0);
scalar surfacez0=0.2;
forAll(mesh.C(),cellI)
{
  // Calculate ustar
  //  vector cc=mesh.C()[cellI];
  scalar x=Foam::pow(1-Foam::neg(MOL)*12*ZREF_/MOL,1.0/3.0);
  scalar psiM=(Foam::pos(MOL)*(-5.0* ZREF_/MOL)+Foam::neg(MOL)*(1.5*Foam::log(1.0/3.0*(1+x+Foam::sqr(x)))-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0))+3.14/Foam::sqrt(3.0)));
  scalar ustar= M*kappa_/(Foam::log((ZREF_+Z0_)/Z0_)-psiM+SMALL);
  // Roughness correction
  ustar=ustar*Foam::log(160.0/Z0_)/Foam::log(160.0/surfacez0);
  // Calculate wind speed 
  scalar z_ = y_[cellI];
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
  U[cellI] = upperstablecorrection*locFlowDir*ustar/kappa_*(Foam::log((z_+surfacez0)/surfacez0)-psiM+extraterm+SMALL);
  // Calculate Turbulence
  scalar phiM=(Foam::pos(MOL)*(1+ 5.0*z_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*12*Foam::pos(z_)*z_/MOL,-1.0/3.0)));
  nut[cellI]=ustar*kappa_*(z_+surfacez0)/phiM*Foam::sqr(1-z_/PBLH);
  alphaT[cellI]=nut[cellI]/Prt.value();
  // Calculate Temperature 
  scalar thetastar=TRef*Foam::sqr(ustar)/(kappa_*9.81*MOL);
  x=Foam::pow(1- Foam::neg(MOL)*16*(TZRef_)/MOL,0.25);
  scalar psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*z_/MOL);
  scalar Twall_= TRef - thetastar/kappa_*(Foam::log((TZRef_+surfacez0)/surfacez0)-psiHL);
  z_=min(z_,0.2*PBLH)+surfacez0;
  x=Foam::pow(1- neg(MOL)*16*(z_)/MOL,0.25);
  psiHL= (MOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x)))): (-5.0*z_/MOL);
  T[cellI]=Twall_+thetastar/kappa_*(Foam::log((z_+ surfacez0)/surfacez0)-psiHL);
}
U.correctBoundaryConditions();
T.correctBoundaryConditions();
U.write();
T.write();
nut.correctBoundaryConditions();
alphaT.correctBoundaryConditions();
nut.write();
alphaT.write();
// Wind
Info<<"Wind Bounds:"<<gMin(U)<<"    "<<gMax(U)<<endl;
// Temperature
Info<<"Potential Temperature Bounds:"<<gMin(T)<<"    "<<gMax(T)<<endl;
// Turbulence
Info<<"Nut:"<<gMin(nut)<<" "<<gMax(nut)<<endl;


