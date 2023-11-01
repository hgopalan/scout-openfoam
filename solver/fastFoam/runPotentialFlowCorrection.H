// Infer the velocity potential BCs from the pressure
wordList PhiBCTypes
(
 p.boundaryField().size(),
 zeroGradientFvPatchScalarField::typeName
 );

forAll(p.boundaryField(), patchi)
{
  if (p.boundaryField()[patchi].fixesValue())
    {
      PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

Info<< "Running Potential Flow Solver\n" << endl;
volScalarField Phi
(
 IOobject
 (
  "Phi",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar(dimLength*dimVelocity, Zero),
 PhiBCTypes
 );

label PhiRefCell = 0;
scalar PhiRefValue = 0;
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
 Phi,
 pimple.dict(),
 pRefCell,
 pRefValue
 );
mesh.setFluxRequired(Phi.name());

phi = fvc::flux(U);

adjustPhi(phi, U, p);

// Non-orthogonal velocity potential corrector loop
scalar counter=0;
while (counter<=5)
  {
    fvScalarMatrix PhiEqn
      (
       fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
       ==
       fvc::div(phi)
       );

    PhiEqn.setReference(PhiRefCell, PhiRefValue);
    PhiEqn.solve();
    phi -= PhiEqn.flux();
    scalar contError=mag(fvc::div(phi))().weightedAverage(mesh.V()).value();
    if (counter>=5 || contError<1e-5)
      {
	break;
      }
    counter=counter+1;
  }
Info<< "Continuity error = "
<< mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
<< endl;

U = fvc::reconstruct(phi);
U.correctBoundaryConditions();

Info<< "Interpolated velocity error = "
<< (sqrt(sum(sqr(fvc::flux(U) - phi)))/sum(mesh.magSf())).value()
<< endl;
Info<<"Wind Bounds:"<<gMin(U)<<"    "<<gMax(U)<<endl;
// Write U
//if(runTime.value()==0)
//  U.write();
Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
<< nl << endl;



