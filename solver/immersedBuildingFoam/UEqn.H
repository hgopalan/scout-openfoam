// Solve the momentum equation
volScalarField sourceU= treeDragCoeff*mag(U);
MRF.correctBoundaryVelocity(U);
// Solve the momentum equation
// Correcting Momentum 
// Damping
volTensorField gradU=fvc::grad(U);
fvVectorMatrix UEqn
(
 fvm::ddt(U) 
 ==
 -fvm::div(phi, U)
 +fvm::laplacian(nut, U)
 +fvc::div(nut*dev2(gradU.T()))
 +nudgeU_
 -fvm::Sp(sourceU,U)
 +fvOptions(U)
 );
UEqn.relax();
fvOptions.constrain(UEqn);
//UEqn.setValues(blankGhostCells,blankGhostWind);
if (pimple.momentumPredictor())
  {
    solve
      (
       UEqn
       ==
       buildingBlanking*fvc::reconstruct
       (
	(
	 - ghf*fvc::snGrad(rhok)
	 - fvc::snGrad(p_rgh)
	 )*mesh.magSf()
	)
       );
    fvOptions.correct(U);
  }



