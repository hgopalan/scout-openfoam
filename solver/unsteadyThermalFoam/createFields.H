Info<<"Reading Fields\n"<<endl;

volScalarField p_rgh
(
 IOobject
 (
  "p_rgh",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

volVectorField U
(
 IOobject
 (
  "U",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

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

volScalarField thetaV
(
 IOobject
 (
  "thetaV",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 T
 );


volScalarField qv
(
 IOobject
 (
  "qv",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );
volScalarField relHum
(
 IOobject
 (
  "relHum",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 0.8+0*mag(U),
 zeroGradientFvPatchField<scalar>::typeName
 );
/*volScalarField presresidual
(
 IOobject
 (
  "presresidual",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 0.0+0*mag(U),
 zeroGradientFvPatchField<scalar>::typeName
 );
volScalarField magU
(
 IOobject
 (
  "magU",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mag(U)
 ); */
// Solver Type 
enum solverType
  {
    steadyWind,
    unsteadyWind,
    steadyThermal,
    unsteadyThermal,
    steadyHumid,
    unsteadyHumid
  };
static const Foam::Enum
<
  solverType
  >
solverTypeNames_
({
  { solverType::steadyWind, "steadyWind" },
    { solverType::steadyThermal, "steadyThermal" },
      { solverType::steadyHumid, "steadyHumid" },
	{ solverType::unsteadyWind, "unsteadyWind" },
	  { solverType::unsteadyThermal, "unsteadyThermal" },
	    { solverType::unsteadyHumid, "unsteadyHumid" },
	      });
// Calculate the Flux 
#include "createPhi.H"
// Read the properties of the atmospheric boundary layer 
#include "readABLProperties.H"
// Setup the trees if need to be included 
#include "readTreeProperties.H"
// Setup Source if need to be included
#include "sourceProperties.H"
// Setup turbulence information 
#include "createTurbulenceFields.H"
// Density modification due to Temperature (Boussinesq correction)
scalar TRef=interpolateXY(U.time().value(),timevalues,TRefvalues);
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
// Gravity related terms (not sure if we need so many?)
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
// Pressure is calculated for post-processing from p_rgh and setting reference value 
volScalarField p
(
 IOobject
 (
  "p",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 p_rgh + rhok*gh
 );
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
 p,
 p_rgh,
 pimple.dict(),
 pRefCell,
 pRefValue
 );
if (p_rgh.needReference())
  {
    p += dimensionedScalar
      (
       "p",
       p.dimensions(),
       pRefValue - getRefCellValue(p, pRefCell)
       );
  }
mesh.setFluxRequired(p_rgh.name());
#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"


volScalarField blankcell
(
 IOobject
 (
  "blankcell",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 0+0.0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );

// Test Code
//#include "identifyUnphysicalCells.H"
#include "initializeSolarCalculator.H"
scalar generateWind=1;
scalar generateThermal=1;


