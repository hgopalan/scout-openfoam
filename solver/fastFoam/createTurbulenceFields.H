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
volScalarField y_ = wallDist(mesh).y();

