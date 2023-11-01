volScalarField Tmrt
(
 IOobject
 (
  "Tmrt",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 T,
 zeroGradientFvPatchField<scalar>::typeName
 );
 volScalarField PMV
(
 IOobject
 (
  "PMV",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 1.0+0.0*mag(U),
 zeroGradientFvPatchField<scalar>::typeName
 );

