Info<< "Reading transportProperties\n" << endl;

IOdictionary rainDict
(
 IOobject
 (
  "rainDict",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
  )
 );

dimensionedScalar Rh("Rh",rainDict);
//dimensionedScalar temp("temp",rainDict);
dimensionedScalar temp=TRef;
dimensionedScalar rhoa
(
 "rhoa",
 dimensionSet(1,-3,0,0,0,0,0),
 scalar(0) 
 );
dimensionedScalar rhop
(
 "rhop",
 dimensionSet(1,-3,0,0,0,0,0),
 scalar(0) 
 );
dimensionedScalar mua
(
 "mua",
 dimensionSet(1,-1,-1,0,0,0,0),
 scalar(0) 
 );
/*scalar scalingFactor
(
 rainDict.lookupOrDefault<scalar>("scalingFactor", 1.0)
 );*/
//dimensionedScalar UMax_("UMax",rainDict);
//dimensionedScalar DIn("DIn",rainDict);
//dimensionedScalar length("length",rainDict);
List<scalar> defaultPhases = List<scalar> (4);
defaultPhases[0]=scalar(0.0005);
defaultPhases[1]=scalar(0.001);
defaultPhases[2]=scalar(0.002);
defaultPhases[3]=scalar(0.005);

List<scalar> phases
(
 rainDict.lookupOrDefault<List<scalar> >("phases",defaultPhases)
 );

List<scalar> fh;
fh.setSize(phases.size());

List<scalar> Vt;
Vt.setSize(phases.size());


volScalarField Re
(
 IOobject
 (
  "Re",
  "0",
  //runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("Re",dimensionSet(0,0,0,0,0,0,0),500)
 );

volScalarField CdRe
(
 IOobject
 (
  "CdRe",
  "0",
  //  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("CdRe",dimensionSet(0,0,0,0,0,0,0),250)
 );
