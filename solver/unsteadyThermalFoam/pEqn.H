{
  volScalarField tempP=p_rgh;
  volScalarField rAU("rAU", 1.0/UEqn.A());
  surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
  volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

  surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());
  surfaceScalarField phiHbyA
    (
     "phiHbyA",
     fvc::flux(HbyA)
     + rAUf*fvc::ddtCorr(U, phi)
     + phig
     );
  MRF.makeRelative(phiHbyA);
  /*if(Steady || multiSteady)
    {
    phiHbyA=fvc::flux(HbyA);
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);
    phiHbyA += phig;
    }*/
  // Update the pressure BCs to ensure flux consistency
  //    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
#include "adjustPressureBC.H"
  // Taken from cgins
  volScalarField divU=0.0/500*fvc::div(U);
  while (pimple.correctNonOrthogonal())
    {
      fvScalarMatrix p_rghEqn
        (
	 fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)+divU
	 );
      p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
      //        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
      //p_rghEqn.diag()+=scalar(1e-15);
      p_rghEqn.solve();
      if (pimple.finalNonOrthogonalIter())
        {
	  // Calculate the conservative fluxes
	  phi = phiHbyA - p_rghEqn.flux();

	  // Explicitly relax pressure for momentum corrector
	  p_rgh.relax();

	  // Correct the momentum source with the pressure gradient flux
	  // calculated from the relaxed pressure
	  U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
	  U.correctBoundaryConditions();
        }
    }

#include "continuityErrs.H"
  // Bad U Check
  forAll(U,cellI)
    {
      if((U[cellI][2]>5 || U[cellI][2]<-5)  && T[cellI]>=308.0)
	{
	  Pout<<"Unrealistic:"<<U[cellI]<<"  "<<T[cellI]<<"   "<<y_[cellI]<<endl;
	  U[cellI][2]=0;
	}
      //      if(mag(U[cellI][0])>5.0 || mag(U[cellI][1])>5.0 || mag(U[cellI][2])>5.0)
      //	Pout<<mesh.C()[cellI]<<endl;
    }
  U.correctBoundaryConditions();
  p = p_rgh + rhok*gh;
  //  magU=mag(U);
  if (p_rgh.needReference())
    {
      p += dimensionedScalar
        (
	 "p",
	 p.dimensions(),
	 pRefValue - getRefCellValue(p, pRefCell)
	 );
      p_rgh = p - rhok*gh;
    }
  //  presresidual=tempP-p_rgh;
  //  presresidual.correctBoundaryConditions();
  //  scalar physicalResidual=0.0;
  //  scalar unphysicalResidual=0.0;
  //  forAll(presresidual,cellI)
  //    {
  //      if(unphysicalCells[cellI]>3)
  //	unphysicalResidual=unphysicalResidual+presresidual[cellI];
  //      else
  //	physicalResidual=physicalResidual+presresidual[cellI];
  //    }
  //  reduce(unphysicalResidual,sumOp<scalar>());
  //  reduce(physicalResidual,sumOp<scalar>());
  treeDragForce=treeDragCoeff*mag(U)*U;
  Info<<"Wind Bounds:"<<gMin(U)<<"    "<<gMax(U)<<endl;      
  //  Info<<"Physical Residual:"<<physicalResidual<<endl;
  //  Info<<"Unphysical Residual:"<<unphysicalResidual<<endl;
  //  Info<<"Normalized Physical Residual:"<<physicalResidual/physicalCounter<<endl;
  //  Info<<"Normalized Unphysical Residual:"<<unphysicalResidual/(unphysicalCounter+SMALL)<<endl;
}
