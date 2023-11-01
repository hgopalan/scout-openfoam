// Turbulent viscosity 
volScalarField nut
(
 IOobject
 (
  "nut",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );
// Thermal Diffusivity 
volScalarField alphaT
(
 IOobject
 (
  "alphaT",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );



dimensionedScalar temp_C1=0.1284;
dimensionedScalar turb_sigma=1.0;
dimensionedScalar turb_kappa=0.41;
dimensionedScalar C2=temp_C1/sqr(turb_kappa)+turb_sigma;
volScalarField y_ = wallDist(mesh).y();
volScalarField turb_C2=C2+0*y_;
Info<<" Minimum Wall Distance:"<<gMin(y_)<<endl;
volScalarField turb_C1
(
 IOobject
 (
  "turb_C1",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 temp_C1+0*nut,
 zeroGradientFvPatchField<scalar>::typeName
 );

volScalarField turb_C3
(
 IOobject
 (
  "turb_C3",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 0*nut,
 zeroGradientFvPatchField<scalar>::typeName
 );

