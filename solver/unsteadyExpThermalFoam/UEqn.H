// Solve the momentum equation
volScalarField sourceU= treeDragCoeff*mag(U);
MRF.correctBoundaryVelocity(U);
// Solve the momentum equation

// Damping
volTensorField gradU=fvc::grad(U);
fvVectorMatrix UEqn
(
 fvm::ddt(U) 
 ==
 -fvm::div(phi, U)
 +fvm::laplacian((nu+nut), U)
 +fvc::div((nu+nut)*dev2(gradU.T()))
 +nudgeU_
 -fvm::Sp(sourceU,U)
 +fvOptions(U)
 + MRF.DDt(U)
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
	 - fvc::snGrad(p_rgh)
	 )*mesh.magSf()
	)
       );

    fvOptions.correct(U);
  }


