{
  volScalarField tempP=p_rgh;
  volScalarField rAU("rAU", 1.0/(1e-5+UEqn.A()));
  surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
  volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
  surfaceScalarField phig(-buildingBlankingSurface*rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());
  surfaceScalarField phiHbyA
    (
     "phiHbyA",
     buildingBlankingSurface*fvc::flux(HbyA)
     + buildingBlankingSurface*rAUf*fvc::ddtCorr(U, phi)
     + phig
     );
  MRF.makeRelative(phiHbyA);
#include "adjustPressureBC.H"
  while (pimple.correctNonOrthogonal())
    {
      fvScalarMatrix p_rghEqn
        (
	 fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
	 );
      p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
      p_rghEqn.solve();
      if (pimple.finalNonOrthogonalIter())
        {
	  // Calculate the conservative fluxes
	  phi = (phiHbyA - p_rghEqn.flux());

	  // Explicitly relax pressure for momentum corrector
	  p_rgh.relax();

	  // Correct the momentum source with the pressure gradient flux
	  // calculated from the relaxed pressure
	  volVectorField field1=fvc::reconstruct(phig);
	  volVectorField field2=fvc::reconstruct(p_rghEqn.flux());
	  U = (HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf));
	  U.correctBoundaryConditions();
	  fvOptions.correct(U);
	  forAll(U,cellI)
	  {
		  if(buildingBlanking[cellI]==0 && mag(U[cellI])>0)
			Pout<<U[cellI]<<"  "<<field1[cellI]<<"  "<<field2[cellI]<<endl;
			
       }	  	  
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
    }
  U.correctBoundaryConditions();
  p = p_rgh + rhok*gh;

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
  presresidual=tempP-p_rgh;
  presresidual.correctBoundaryConditions();
  scalar physicalResidual=0.0;
  scalar unphysicalResidual=0.0;
  forAll(presresidual,cellI)
    {
      if(unphysicalCells[cellI]>3)
	unphysicalResidual=unphysicalResidual+presresidual[cellI];
      else
	physicalResidual=physicalResidual+presresidual[cellI];
    }
  reduce(unphysicalResidual,sumOp<scalar>());
  reduce(physicalResidual,sumOp<scalar>());
  Info<<"Wind Bounds:"<<gMin(U)<<"    "<<gMax(U)<<endl;      
  Info<<"Physical Residual:"<<physicalResidual<<endl;
  Info<<"Unphysical Residual:"<<unphysicalResidual<<endl;
  Info<<"Normalized Physical Residual:"<<physicalResidual/physicalCounter<<endl;
  Info<<"Normalized Unphysical Residual:"<<unphysicalResidual/(unphysicalCounter+SMALL)<<endl;
}

