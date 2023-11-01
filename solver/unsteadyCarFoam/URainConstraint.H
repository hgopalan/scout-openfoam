//  const scalarField& vol = mesh.V();
label nDamped = 0;
//  volScalarField umag = mag((Urain[phase_no]));
volScalarField umag = mag((UCheat));
scalar umax = max(umag).value();
umax=Foam::sqrt(Foam::pow(umax,2)+Foam::pow(max(Vt),2));
forAll(Urain[phase_no], cellI)
{
  scalar magU = mag(Urain[phase_no][cellI]);
  //      if (magU > UMax_.value())
  if (magU > 2.5*umax)
    {
      //	  scalar scale = sqr(Foam::cbrt(vol[cellI]));
      //	  UrainEqn.diag()[cellI] += scale*(magU-UMax_.value());
      //	  UrainEqn.diag()[cellI] += scale*(magU-umax);
      Urain[phase_no][cellI] = ((2.5*umax)/magU)*Urain[phase_no][cellI];
      nDamped++;
    }
}
reduce(nDamped, sumOp<label>());
Info<< nDamped << " ("
<< 100*scalar(nDamped)/mesh.globalData().nTotalCells()
<< "%) of cells" << endl;

