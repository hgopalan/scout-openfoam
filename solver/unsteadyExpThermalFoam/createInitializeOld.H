Info<<"TRef:"<<TRef<<endl;
forAll(mesh.C(),i)
{
  scalar kappa_ = 0.41;
  scalar locMOL=interpolateXY(U.time().value(),timevalues,movalues);
  scalar locU1=interpolateXY(U.time().value(),timevalues,uvalues);
  scalar locV1=interpolateXY(U.time().value(),timevalues,vvalues);
  scalar locM=Foam::sqrt(sqr(locU1)+sqr(locV1));  
  scalar z_ = y_[i];
  z_=min(z_,0.9*PBLH.value())+Z0_;
  //  scalar xZ_ = Foam::pow(1- neg(locMOL)*16*ZREF_/locMOL,0.25);
  scalar x=Foam::pow(1-Foam::neg(MOL.value())*12*ZREF_/MOL.value(),1.0/3.0);
  //  scalar intermediateZ_ = (locMOL < 0) ? (Foam::log(0.5*(1+sqr(xZ_))*sqr(0.5*(1+xZ_)))-2*Foam::atan(xZ_)+0.5*3.14) : (-5.0* ZREF_/locMOL) ;
  scalar intermediateZ_ = (locMOL < 0) ? (1.5*Foam::log(1.0/3.0*(1+x+Foam::sqr(x)))-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0))+3.14/Foam::sqrt(3.0)) : (-5.0* ZREF_/locMOL) ;
  //  scalar x_ = Foam::pow(1- neg(locMOL)*16*z_/locMOL,0.25);
  //  scalar intermediateU_ = (locMOL < 0) ? (Foam::log(0.5*(1+sqr(x_))*sqr(0.5*(1+x_)))-2*Foam::atan(x_)+0.5*3.14) : (-5.0* z_/locMOL) ;
  x=Foam::pow(1-Foam::neg(MOL.value())*12*z_/MOL.value(),1.0/3.0);
  scalar intermediateU_ = (locMOL < 0) ? (1.5*Foam::log(1.0/3.0*(1+x+Foam::sqr(x)))-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0))+3.14/Foam::sqrt(3.0)) : (-5.0* z_/locMOL) ;
  vector locFlowDir(locU1/locM,locV1/locM,0);
  scalar locUstar=locM*kappa_/(Foam::log((ZREF_+Z0_)/Z0_)-intermediateZ_+SMALL); 
  scalar extraterm=0;
  scalar upperstablecorrection=1;
  if(mag(MOL.value())>500)
    extraterm=z_/MBL.value()-z_/PBLH.value()*0.5*z_/MBL.value();    
  else if(pos(MOL.value()))
    {
      if(y_[i]>PBLH.value())
	upperstablecorrection=Foam::min(y_[i]/(0.9*PBLH.value()),1.2);
      extraterm=z_/MBL.value()-z_/PBLH.value()*0.5*z_/MBL.value()-5*z_/MOL.value()*0.5*z_/PBLH.value();    
    }
  else
    {
      //      if(y_[i]>PBLH.value())
      //        upperstablecorrection=Foam::min(y_[i]/(0.9*PBLH.value()),1.2);
      extraterm=z_/MBL.value()-z_/PBLH.value()*0.5*z_/MBL.value()+z_/PBLH.value()*(1+(Foam::pow(1-12*z_/MOL.value(),2.0/3.0)-1)/(8*z_/MOL.value()));
    }
  U[i] = upperstablecorrection*locFlowDir*locUstar/kappa_*(Foam::log((z_+Z0_)/Z0_)-intermediateU_+extraterm+SMALL);
  //  Pout<<z_<<"   "<<locUstar/kappa_*(Foam::log((z_+Z0_)/Z0_)-intermediateU_)<<"  "<<locUstar/kappa_*term2<<"   "<<locUstar/kappa_*term3<<endl;
  scalar thetastar=TRef.value()*Foam::sqr(locUstar)/(kappa_*9.81*locMOL);
  scalar x_=Foam::pow(1- neg(locMOL)*16*(TZRef_)/locMOL,0.25);
  scalar psiHL= (locMOL<0)? (2*Foam::log(0.5*(1+Foam::sqr(x_)))): (-5.0*z_/locMOL);
  scalar Twall_= TRef.value() - thetastar/kappa_*(Foam::log((TZRef_+Z0_)/Z0_)-psiHL);
  x_=Foam::pow(1- neg(locMOL)*16*(z_)/locMOL,0.25);
  scalar intermediateNut_=(locMOL<0) ? Foam::pow(1-16*z_/locMOL,-0.25) : (1+5.0*z_/locMOL) ;
  nut[i]=locUstar*kappa_*(z_+Z0_)/intermediateNut_*Foam::sqr(1-z_/PBLH.value());
  alphaT[i]=nut[i]/Prt.value();
  z_=min(z_,0.2*PBLH.value())+Z0_;
  T[i]=Twall_+thetastar/kappa_*(Foam::log((z_+ Z0_)/Z0_)-psiHL);//+1e-8*pos(y_[i]-PBLH.value())*(y_[i]-PBLH.value());
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


