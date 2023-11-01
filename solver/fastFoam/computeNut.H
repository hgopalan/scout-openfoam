scalar surfaceZ0=0.21;
scalar kappa_=0.41;
forAll(mesh.C(),cellI)
{
  // Calculate ustar
  //  vector cc=mesh.C()[cellI];
  scalar M=mag(U[cellI]);
  scalar z_ = y_[cellI];
  z_=min(z_,0.9*PBLH)+surfaceZ0;
  scalar   x=Foam::pow(1-Foam::neg(MOL)*12*z_/MOL,1.0/3.0);
  scalar psiM=(Foam::pos(MOL)*(-5.0* z_/MOL)+Foam::neg(MOL)*(1.5*Foam::log(1.0/3.0*(1+x+Foam::sqr(x)))-Foam::sqrt(3.0)*Foam::atan((1+2*x)/Foam::sqrt(3.0))+3.14/Foam::sqrt(3.0)));
  scalar ustar= M*kappa_/(Foam::log((z_+surfaceZ0)/surfaceZ0)-psiM+SMALL);
  scalar phiM=(Foam::pos(MOL)*(1+ 5.0*z_/MOL)+Foam::neg(MOL)*(Foam::pow(1- Foam::neg(MOL)*12*Foam::pos(z_)*z_/MOL,-1.0/3.0)));
  nut[cellI]=ustar*kappa_*(z_+surfaceZ0)/phiM*Foam::sqr(1-z_/PBLH);
  alphaT[cellI]=nut[cellI]/Prt.value();
}
nut.correctBoundaryConditions();
alphaT.correctBoundaryConditions();
// Turbulence
Info<<"Nut:"<<gMin(nut)<<" "<<gMax(nut)<<endl;


