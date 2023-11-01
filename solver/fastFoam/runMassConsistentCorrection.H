{
  Info<<" Running Fast Solver Correction \n "<<endl;
  forAll(alpha,cellI)
    {
      alpha[cellI][0]=1.0;
      alpha[cellI][4]=1.0;
      alpha[cellI][8]=alpha[cellI][0]/Foam::pow(mesh.V()[cellI],1.0/10.0);
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

  //  U.write();
  Info<< "Interpolated velocity error = "
      << (sqrt(sum(sqr(fvc::flux(U) - phi)))/sum(mesh.magSf())).value()
      << endl;
  Info<<"Wind Bounds:"<<gMin(U)<<"    "<<gMax(U)<<endl;
  phi=fvc::flux(U);
  Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s"
      << nl << endl;
}
