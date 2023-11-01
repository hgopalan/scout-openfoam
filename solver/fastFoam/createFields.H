Info<< "Reading thermophysical properties\n" << endl;
Info<< "Reading field U\n" << endl;
volVectorField U
(
 IOobject
 (
  "U",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );
Info<< "Reading field T\n" << endl;
volScalarField T
(
 IOobject
 (
  "T",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );
volScalarField p
(
 IOobject
 (
  "p",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );
volScalarField Rib
(
 IOobject
 (
  "Rib",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 1.225+0.0*mag(U),
 zeroGradientFvPatchField<scalar>::typeName
 );

 
#include "createPhi.H"

#include "readABLProperties.H"
#include "readTreeProperties.H"
#include "createTurbulenceFields.H"
//Info<< "Creating turbulence model\n" << endl;
//autoPtr<incompressible::turbulenceModel> turbulence
//(
// incompressible::turbulenceModel::New(U, phi, laminarTransport)
// );
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
// Kinematic density for buoyancy force
volScalarField rhok
(
 IOobject
 (
  "rhok",
  runTime.timeName(),
  mesh
  ),
 1.0 - 1/TRef*(T - TRef)
 );

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField buildingZone
(
 IOobject
 (
  "buildingZone",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  ),
 mesh
 );


#include "createIncompressibleRadiationModel.H"
volScalarField blankcell
(
 IOobject
 (
  "blankcell",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 0.0*mag(U),
 zeroGradientFvPatchField<scalar>::typeName
 );
if(! multiSteady)
  {
    Info<<"Initializing...."<<endl;
#include "createInitialize.H"
    //#include "runPotentialFlowCorrection.H"
#include "wallBlocking.H"
  }


// Sponge Layer 
volScalarField spongeLayerViscosity
(
 IOobject
 (
  "spongeLayerViscosity",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("spongeLayerViscosity",dimensionSet(0,1, -1, 0, 0, 0, 0),0.0),
 zeroGradientFvPatchField<scalar>::typeName
 );

scalar ztop=1200;
scalar zbottom=800;
forAll(U,cellI)
{
  if(y_[cellI]>zbottom)
    spongeLayerViscosity[cellI]=Foam::sqr(Foam::sin(3.14/2*(1-(ztop-y_[cellI])/(ztop-zbottom))));
}
spongeLayerViscosity.correctBoundaryConditions();
spongeLayerViscosity.write();

volScalarField lambda
(
 IOobject
 (
  "lambda",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  ),
 mesh
 );
volTensorField alpha
(
 IOobject
 (
  "alpha",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 0*fvc::grad(U),
 zeroGradientFvPatchField<tensor>::typeName
 );
#include "readScalarProperties.H"
#include "createScalarFields.H"
