// Neutral 
Info<<"PBLH:"<<PBLH<<"   "<<MBL<<endl;
// Not sure if wall distance or cell center should be used?
volScalarField z_=min(y_,0*y_+0.9*PBLH);//+Z0_;
//volScalarField z_=min((mesh.C()).component(vector::Z),0.9*PBLH+0*y_)+Z0_;
tmp<volScalarField> R=turb_kappa*z_*Foam::sqr(1-z_/PBLH);
tmp<volScalarField> S=1/turb_kappa*(1-z_/PBLH)*(1/z_+1/MBL+1/(PBLH-z_));
tmp<volScalarField> prod=R*S;
tmp<volScalarField> R1=turb_kappa*z_*Foam::sqr(1-z_/PBLH);
tmp<volScalarField> S1=1/turb_kappa*(1-z_/PBLH)*(1/z_+1/MBL+1/(PBLH-z_));
tmp<volScalarField> dSdz=-1/(turb_kappa*PBLH)*(1/z_+1/MBL+1/(PBLH-z_))+1/turb_kappa*(1-z_/PBLH)*(-1/(z_*z_));
tmp<volScalarField> diss=turb_C2*Foam::sqr(R1/S1*dSdz);
tmp<volScalarField> dRdz=turb_kappa*(1+3*Foam::sqr(z_/PBLH)-4*z_/PBLH);
tmp<volScalarField> dR2dz2=turb_kappa*(6*z_/Foam::sqr(PBLH)-4/PBLH);
tmp<volScalarField> diff=turb_sigma*(Foam::sqr(dRdz)+turb_kappa*z_*Foam::sqr(1-z_/PBLH)*dR2dz2);
turb_C1=(diss-diff)/prod;
turb_C1.correctBoundaryConditions();
Info<<"C1:"<<max(turb_C1).value()<<"  "<<min(turb_C1).value()<<endl;
if(Foam::mag(MOL)>500)
  {
    turb_C3=0*turb_C3;
    Info<<"C3:"<<max(turb_C3).value()<<"  "<<min(turb_C3).value()<<endl;
  }
 else if(pos(MOL) && MOL<2.5)
   {
     volScalarField R=turb_kappa*z_/6.0*Foam::sqr(1-z_/PBLH);
     volScalarField S=(1-z_/PBLH)/turb_kappa*(6/z_+1/MBL+1/(PBLH-z_));
     tmp<volScalarField> prod=turb_C1*R*S;
     tmp<volScalarField> dSdz=1/turb_kappa*(-1/PBLH*(6/z_+1/MBL+1/(PBLH-z_))+(1-z_/PBLH)*(-6/(z_*z_)));
     tmp<volScalarField> diss=turb_C2*Foam::sqr(R/S*dSdz);
     tmp<volScalarField> dRdz=turb_kappa/6.0*(1-4*z_/PBLH+3*Foam::sqr(z_/PBLH));
     tmp<volScalarField> dR2dz2=turb_kappa/6*(6*z_/(PBLH*PBLH)-4/PBLH);
     tmp<volScalarField> diff=turb_sigma*(Foam::sqr(dRdz)+R*dR2dz2);
     turb_C3=1/(R*S)*(prod-diss+diff);
     turb_C3.correctBoundaryConditions();
     Info<<"C3:"<<max(turb_C3).value()<<"  "<<min(turb_C3).value()<<endl;
   }
 else if(pos(MOL))
   {
     volScalarField phiM=1+5*z_/MOL;
     volScalarField fmu=(1-z_/PBLH);
     volScalarField alpha1=phiM+z_/MBL*fmu+z_/PBLH*(1-phiM);
     volScalarField alpha2=alpha1-(5*z_/MOL+z_/MBL*fmu-Foam::sqr(z_)/(MBL*PBLH)-10*Foam::sqr(z_)/(MOL*PBLH));
     tmp<volScalarField> prod=turb_C1*Foam::sqr(fmu)/phiM*alpha1;     
     tmp<volScalarField> diss=turb_C2*sqr(turb_kappa/phiM*alpha2/alpha1*Foam::sqr(fmu));
     tmp<volScalarField> R=turb_kappa*z_/phiM*Foam::sqr(fmu);
     tmp<volScalarField> dRdz=turb_kappa*fmu/Foam::sqr(phiM)*(1-3*z_/PBLH-10*Foam::sqr(z_)/(PBLH*MOL));
     tmp<volScalarField> dR2dz2=turb_kappa/Foam::pow(phiM,3)*(6*z_/Foam::sqr(PBLH)+30/MOL*Foam::sqr(z_/PBLH)+50*Foam::pow(z_,3)/Foam::sqr(PBLH*MOL)-4/PBLH-10/MOL);
     tmp<volScalarField> diff=turb_sigma*(Foam::sqr(dRdz)+R*dR2dz2);
     turb_C3=phiM/(alpha1*Foam::sqr(fmu))*(prod-diss+diff);
     turb_C3.correctBoundaryConditions();
     Info<<"C3:"<<max(turb_C3).value()<<"  "<<min(turb_C3).value()<<endl;
   }
 else
   {
	 volScalarField fmu=(1-z_/PBLH);
     volScalarField phiM=Foam::pow(1-12*z_/MOL,-1.0/3.0);
     volScalarField alpha1=phiM+z_/MBL*fmu+z_/PBLH*(1-phiM);
     volScalarField alpha2=alpha1-(4*z_/MOL*Foam::pow(phiM,4)+z_/MBL*(1-z_/PBLH)-Foam::sqr(z_)/(MBL*PBLH)+z_/PBLH*(1-phiM)-4*Foam::sqr(z_)/(PBLH*MOL)*Foam::pow(phiM,4));
     tmp<volScalarField> prod=turb_C1*Foam::sqr(fmu)/phiM*alpha1;
     tmp<volScalarField> diss=turb_C2*sqr(turb_kappa/phiM*alpha2/alpha1*Foam::sqr(fmu));
     tmp<volScalarField> R=turb_kappa*z_/phiM*Foam::sqr(fmu);
     tmp<volScalarField> dRdz=turb_kappa*(Foam::sqr(fmu)/phiM-4*z_/MOL*Foam::sqr(fmu*phiM)-2*z_/PBLH*fmu/phiM);
     tmp<volScalarField> dR2dz2=turb_kappa*(2*z_/Foam::sqr(PBLH)*1/phiM-32*z_*Foam::sqr(fmu/MOL)*Foam::pow(phiM,5)-8/MOL*Foam::sqr(fmu*phiM)-4*fmu/PBLH*1/phiM+16*z_/(MOL*PBLH)*fmu*Foam::sqr(phiM));
     tmp<volScalarField> diff=turb_sigma*(Foam::sqr(dRdz)+R*dR2dz2);
     turb_C3=phiM/(alpha1*Foam::sqr(fmu))*(prod-diss+diff);
     turb_C3.correctBoundaryConditions();
     Info<<"C3:"<<max(turb_C3).value()<<"  "<<min(turb_C3).value()<<endl;
   }




/*{
// Calculate C1 
volScalarField z_=min(y_,0*y_+0.9*PBLH)+Z0_;
volScalarField fmu=sqr(1-z_/PBLH);
//turb_C2=1/fmu*(turb_C1/sqr(turb_kappa)+turb_sigma*(fmu+10*Foam::sqr(z_/PBLH)-4*z_/PBLH*(1+Foam::sqrt(fmu))));
volScalarField alpha1=1+z_/MBL*(1-z_/PBLH);
volScalarField alpha2=1+Foam::sqr(z_)/(PBLH*MBL);
turb_C1=Foam::sqr(turb_kappa)/alpha1*(turb_C2*fmu*Foam::sqr(alpha2/alpha1)-turb_sigma*(fmu+10*Foam::sqr(z_/PBLH)-4*z_/PBLH*(1+Foam::sqrt(fmu))));
if(neg(MOL.value()))
{
volScalarField phiM=Foam::pow(1-16*z_/MOL,-0.25);
tmp<volScalarField> prod=turb_C1*fmu;
tmp<volScalarField> diss=turb_C2*Foam::sqr(turb_kappa*fmu/phiM*(1-4*z_/MOL*Foam::pow(phiM,4)));
volScalarField R=turb_kappa*z_/phiM*fmu;
volScalarField dRdz=turb_kappa*(fmu/phiM-4*z_/MOL*Foam::pow(phiM,3)-2*z_/PBLH*Foam::sqrt(fmu)*1/phiM);
volScalarField dR2dz2=turb_kappa*(2*z_/Foam::sqr(PBLH)*1/phiM-48*z_/Foam::sqr(MOL)*fmu*Foam::pow(phiM,7));
dR2dz2=dR2dz2-turb_kappa*(8/MOL*fmu*Foam::pow(phiM,3)+4/PBLH*(1-z_/PBLH)*1/phiM-16*z_/PBLH*1/MOL*Foam::sqrt(fmu)*Foam::pow(phiM,3));
tmp<volScalarField> diff=turb_sigma*(Foam::sqr(dRdz+R*dR2dz2));
turb_C3=(prod-diss+diff);
}
else
{
volScalarField phiM=1+5*z_/MOL;
volScalarField alpha3=(1+z_/PBLH*(1/phiM-1)+z_/(MBL*phiM)*(1-z_/PBLH));
volScalarField alpha4=1+5*Foam::sqr(z_)/PBLH*(5/MOL+1/MBL);
tmp<volScalarField> prod=turb_C1+0*nut;
tmp<volScalarField> diss=turb_C2*fmu*Foam::sqr(turb_kappa)/Foam::pow(phiM,4);
volScalarField R_=turb_kappa*z_/phiM;
volScalarField dRdz=turb_sigma*Foam::sqr(turb_kappa/phiM)*(1-3*z_/PBLH-10*Foam::sqr(z_)/(PBLH*MOL));
volScalarField dR2dz2=2.0/5.0*MOL/Foam::sqr(PBLH)-2/(5*MOL*Foam::pow(phiM,3))*(25+10*MOL/PBLH+Foam::sqr(MOL/PBLH));
dR2dz2=turb_sigma*turb_kappa*dR2dz2;
tmp<volScalarField> diffusion=Foam::sqr(dRdz)+R_*dR2dz2;
turb_C3=Foam::neg(MOL.value()-500)/alpha3*(prod*alpha3-diss*Foam::sqr(alpha4/alpha3)+diffusion);
} 
turb_C3.correctBoundaryConditions();
}*/
