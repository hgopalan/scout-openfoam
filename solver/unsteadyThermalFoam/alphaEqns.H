
fvScalarMatrix alphaEqn
(
 fvm::div(phirain[phase_no],alpharain[phase_no],"div(phirain,alpharain)")
);

alphaEqn.relax();
forAll(alpharain[phase_no], cellI) {
      //alpharain[phase_no][cellI]+=scalar(1e-15);
      alphaEqn.diag()[cellI] += scalar(1e-15);
}
alphaEqn.solve();
bound(alpharain[phase_no],SMALL);
forAll(alpharain[phase_no], cellI) 
{
  if(treeblanking[cellI]>0)
    alpharain[phase_no][cellI]=1e-15;
}

alpharain[phase_no].correctBoundaryConditions();


