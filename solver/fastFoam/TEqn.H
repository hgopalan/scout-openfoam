{
  dimensionedScalar Tmin_=283;
  alphaT = nut/Prt;
  alphaT.correctBoundaryConditions();
  volScalarField alphaEff("alphaEff", alphaT);
  // Convert to Potential Temperature
  fvScalarMatrix TEqn
    (
     fvm::ddt(T)
     + fvm::div(phi, T)

     ==
     fvm::laplacian(alphaEff, T)
     +treeLAD*radiation->ST(rhok*cpRef, T)
     );
  TEqn.relax();
  TEqn.solve();
  // In the unphysical gaps in buildings temperature can built-up.
  // Limiting avoids the built up. 
  forAll(T,cellI)
    {
      if(T[cellI]<283.0)
	{
	  counter=counter+1;
	  T[cellI]=TRef;
	}
      if(T[cellI]>325)
	{
	  T[cellI]=TRef;
	  counter=counter+1;
	}
    }
  //  bound(T,Tmin_);
  T.correctBoundaryConditions();
  reduce(counter, sumOp<scalar>());
  //  if(pimple.finalIter())                                                                                                             
  //    radiation->correct(); 
  //  radiation->correct();
  rhok= 1 - 1.0/TRef*(T-TRef);
}
 
