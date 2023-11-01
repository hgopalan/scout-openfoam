{
  Info<<" Running Fast Solver Correction \n "<<endl;
  wordList PhiBCTypes
    (
     p.boundaryField().size(),
     zeroGradientFvPatchScalarField::typeName
     );
  //  scalar locU1=interpolateXY(U.time().value(),timevalues,uvalues);
  //  scalar locV1=interpolateXY(U.time().value(),timevalues,vvalues);

  forAll(p.boundaryField(), patchi)
    {
      // if (p.boundaryField()[patchi].fixesValue())
      if(p.boundaryField()[patchi].patch().name()=="west")
	{
	  PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
	}

      if(p.boundaryField()[patchi].patch().name()=="east")
	{
	  PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
	}

      if(p.boundaryField()[patchi].patch().name()=="south")
	{
	  PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
	}
 
      if(p.boundaryField()[patchi].patch().name()=="north")
	{
	  PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
	}
 
 
      if(p.boundaryField()[patchi].patch().name()=="upper")
	{
	  PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
	}
 
    }

 
  volScalarField lambda
    (
     IOobject
     (
      "lambda",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedScalar(dimLength*dimVelocity, Zero),
     PhiBCTypes
     );
  volTensorField alpha
    (
     IOobject
     (
      "lambda",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     0*fvc::grad(U),
     zeroGradientFvPatchField<tensor>::typeName
     );

  forAll(alpha,cellI)
    {
      alpha[cellI][0]=1.0;
      alpha[cellI][4]=1.0;
      if(y_[cellI]<=5)
	alpha[cellI][8]=1000;
      else
	alpha[cellI][8]=100.0;
    }
  alpha.correctBoundaryConditions();
  fvScalarMatrix lambdaEqn
    (
     fvm::laplacian(alpha,lambda) ==-2*alpha.component(tensor::XX)*fvc::div(U)
     );
  lambdaEqn.relax();
  lambdaEqn.solve();
  volVectorField glambda=fvc::grad(lambda);
  forAll(U,celli)
    {
      U[celli][0]=U[celli][0]+0.5/alpha[celli][0]*glambda[celli][0];
      U[celli][1]=U[celli][1]+0.5/alpha[celli][4]*glambda[celli][1];
      U[celli][2]=U[celli][2]+0.5/alpha[celli][8]*glambda[celli][2];
    }
  U.correctBoundaryConditions();
  Info<< "Interpolated velocity error = "
      << (sqrt(sum(sqr(fvc::flux(U) - phi)))/sum(mesh.magSf())).value()
      << endl;
  Info<<"Wind Bounds:"<<gMin(U)<<"    "<<gMax(U)<<endl;
  Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s"
      << nl << endl;
}
