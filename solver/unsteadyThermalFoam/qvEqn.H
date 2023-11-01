{
  alphaT = nut/Prt+nu;
  alphaT.correctBoundaryConditions();
  volScalarField alphaEff("alphaEff", alphaT);
  // Convert to Potential Temperature
  fvScalarMatrix qvEqn
    (
     fvm::ddt(qv)
     + fvm::div(phi, qv)

     ==
     fvm::laplacian(alphaEff, qv)
     +treeLAD*leafTemperature
     );
  qvEqn.relax();
  qvEqn.solve();
  bound(qv,qVUpper);
  qv.correctBoundaryConditions();
  volScalarField esat=610.094*Foam::exp(17.625*(T-273.15)/(T-30.11));
  volScalarField qsat=(0.622 *esat/(101325-0.378*esat));
  relHum=min(qv/qsat,1.0);
  relHum.correctBoundaryConditions();
  qv=relHum*qsat;
  qv.correctBoundaryConditions();
}
Info<<"Specific Humidity Bounds:"<<gMin(qv)<<"  "<<gMax(qv)<<endl;                                                 
Info<<"Relative  Humidity Bounds:"<<gMin(relHum)<<"  "<<gMax(relHum)<<endl;   
