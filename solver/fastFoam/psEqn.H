if(runScalar)
 {
    dimensionedScalar sc_=1.0;
    dimensionedScalar scT_=0.7;
    for (int i = 0; i < scalarTable.size(); i++)
      {
	volScalarField alphaEff("alphaEff", nut/scT_);
	fvScalarMatrix psEqn
	  (
	   fvm::ddt(ps[i])
	   + fvm::div(phi,ps[i])
	   - fvm::laplacian(alphaEff,ps[i])
	   ==
	   Foam::pos(finalRelTime_[i]-mesh.time().value())*fvOptions(ps[i])
	   );
	psEqn.relax();
	fvOptions.constrain(psEqn);
	psEqn.solve();
	fvOptions.correct(ps[i]);
	ps[i]=pos(ps[i])*ps[i];	
	ps[i].correctBoundaryConditions();
	Info<<pos(finalRelTime_[i]-mesh.time().value())<<"  Scalar Bounds:"<<gMin(ps[i])<<"    "<<gMax(ps[i])<<endl;
      }
  }
