// Solve the momentum equation

volScalarField sourceU= treeDragCoeff*mag(U);
// Solve the momentum equation
volScalarField uz=U.component(vector::Z);
volScalarField fSponge= spongeLayerViscosity*uz*0.2;
volTensorField gradU=fvc::grad(U);
vector zaxis(0,0,1);
fvVectorMatrix UEqn
(
 fvm::ddt(U) + fvm::div(phi, U)
 ==
 fvm::laplacian(nut, U)
 +fvc::div(nut*dev2(gradU.T()))
 -zaxis*fSponge
 -fvm::Sp(sourceU,U)
 +fvOptions(U)
 );
// - fvc::div(nut*dev2(gradU.T())
 //- fvc::div(nuEff_*dev(T(fvc::grad(U))))
 //
UEqn.relax();
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
  {
    solve
      (
       UEqn
       ==
       fvc::reconstruct
       (
	(
	 - ghf*fvc::snGrad(rhok)
	 )*mesh.magSf()
	)
       );

    fvOptions.correct(U);
  }


