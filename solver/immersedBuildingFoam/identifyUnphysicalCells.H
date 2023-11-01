volScalarField unphysicalCells
(
 IOobject
 (
  "unphysicalCells",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 -100+0*mag(U),
 zeroGradientFvPatchField<scalar>::typeName
 );
scalar unphysicalCounter=0;
scalar physicalCounter=0;
const polyBoundaryMesh& patches = mesh.boundaryMesh();
forAll(U,cellI)
{
  const labelList& faces=mesh.cells()[cellI];
  scalar counter=0;
  forAll(faces,faceI)
    {
      label patchId=mesh.boundaryMesh().whichPatch(faces[faceI]);
      if(patchId>-1 && isA<wallPolyPatch>(patches[patchId]))
	counter=counter+1;
    }
  unphysicalCells[cellI]=counter;
  if(counter>3)
    unphysicalCounter=unphysicalCounter+1;
  else
    physicalCounter=physicalCounter+1;
    //	Pout<<mesh.C()[cellI]<<"   "<<counter<<endl;
}
unphysicalCells.write();
reduce(unphysicalCounter,sumOp<scalar>());
reduce(physicalCounter,sumOp<scalar>());
Info<<"Unphysical Cells:"<<unphysicalCounter<<endl;
