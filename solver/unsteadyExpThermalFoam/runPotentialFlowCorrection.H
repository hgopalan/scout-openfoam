// Infer the velocity potential BCs from the pressure
wordList PhiBCTypes
(
 p.boundaryField().size(),
 zeroGradientFvPatchScalarField::typeName
 );
scalar locU1=interpolateXY(U.time().value(),timevalues,uvalues);
scalar locV1=interpolateXY(U.time().value(),timevalues,vvalues);

forAll(p.boundaryField(), patchi)
{
  // if (p.boundaryField()[patchi].fixesValue())
  if(p.boundaryField()[patchi].patch().name()=="west" && locU1>0)
    {
      PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }

  if(p.boundaryField()[patchi].patch().name()=="east" && locU1<0)
    {
      PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }

  if(p.boundaryField()[patchi].patch().name()=="south" && locV1>0)
    {
      PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
 
  if(p.boundaryField()[patchi].patch().name()=="north" && locV1<0)
    {
      PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
 
 
  if(p.boundaryField()[patchi].patch().name()=="upper")
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
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar(dimLength*dimVelocity, Zero),
 PhiBCTypes
 );

label PhiRefCell = 0;
scalar PhiRefValue = 0;
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
scalar potentialCounter=0;
// Non-orthogonal velocity potential corrector loop
while (potentialCounter<=10)
  {
    fvScalarMatrix PhiEqn
      (
       fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
       ==
       fvc::div(phi)
       );

    PhiEqn.setReference(PhiRefCell, PhiRefValue);
    PhiEqn.solve();

    if (potentialCounter==10)
      {
	phi -= PhiEqn.flux();
      }
    potentialCounter=potentialCounter+1;
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
if(runTime.value()==0)
  U.write();
Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
<< nl << endl;



