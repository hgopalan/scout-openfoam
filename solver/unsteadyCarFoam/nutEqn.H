volScalarField S_=sqrt(2*(symm(fvc::grad(U)) && fvc::grad(U)))+1e-15;
tmp<volScalarField> gradTerm=(U & fvc::grad(Foam::mag(U)));
tmp<volScalarField> buoyancyCoeff_=-turb_C3*S_;
tmp<volScalarField> cs_= (fvc::grad(S_) & fvc::grad(S_));
tmp<volScalarField> treeTerm_=treeDragCoeff*2.0*Foam::sqrt(0.09)/S_*gradTerm;
volScalarField PTraffic=0*nut;
#include "carTurbulence.H"
tmp<fvScalarMatrix> nutEqn
(
 fvm::ddt(nut)
 + fvm::div(phi, nut)
 ==
 fvm::Sp(turb_C1*S_,nut)
 +fvm::SuSp(buoyancyCoeff_,nut)-fvm::Sp(turb_C2*nut/Foam::sqr(S_)*cs_,nut)
 +fvm::laplacian(turb_sigma*nut+nu,nut)
 -fvm::Sp(treeDragCoeff*8.0/3.0*Foam::mag(U),nut)
 +nudgeNut_
 +PTraffic
 );
nutEqn.ref().relax();
solve(nutEqn);
// +fvm::SuSp(treeTerm_,nut)
//bound(nut,1e-3);
//nut.correctBoundaryConditions();
// Compute Turbulence
scalar nutKaput=0;
scalar zcutoff=0.985*PBLH;
scalar nutcutoff=ustar*kappa_*zcutoff*Foam::sqr(1-zcutoff/PBLH);
nutcutoff=1e-10;
if (pos(MOL) && MOL<75)
  nutcutoff=nutcutoff/6.0;
else if (pos(MOL))
  nutcutoff=0.0;
else if (neg(MOL))
  nutcutoff=nutcutoff;
Info<<"Cut off:"<<nutcutoff<<endl;
forAll(nut,cellI)
{
  vector cc=mesh.C()[cellI];
  scalar z_ = cc[2];
  z_=min(z_,0.9*PBLH);
  scalar correction=Foam::sqr(1-z_/PBLH);
  if(mag(MOL)>500 && nut[cellI]<nutcutoff)
    {
      nut[cellI]=ustar*kappa_*z_*correction;
      nutKaput=nutKaput+1;
    }
  else if (pos(MOL) && MOL<75 && nut[cellI]<nutcutoff)
    {
      nut[cellI]=ustar*kappa_*z_/6*correction;
      nutKaput=nutKaput+1;
    }
  else if (pos(MOL) && nut[cellI]<nutcutoff)
    {
      nut[cellI]=ustar*kappa_*z_/(1+5*z_/MOL)*correction;
      nutKaput=nutKaput+1;
    }
   else if (neg(MOL) && nut[cellI]<nutcutoff)
   {
     nut[cellI]=ustar*kappa_*z_/Foam::pow(1- 12*z_/MOL,-1.0/3.0)*correction;
     //     nut[cellI]=nutcutoff;
     nutKaput=nutKaput+1;
   }
}
reduce(nutKaput, sumOp<scalar>());
Info<<"Nut:"<<gMin(nut)<<" "<<gMax(nut)<<endl;   
Info<<"Out of bound nut cells:"<<nutKaput<<endl;
