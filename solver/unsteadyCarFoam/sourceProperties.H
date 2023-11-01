// Velocity Source 
volVectorField velocitySource
(
 IOobject
 (
  "velocitySource",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 0*U,
 zeroGradientFvPatchField<scalar>::typeName
 );
volScalarField temperatureSource
(
 IOobject
 (
  "temperatureSource",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
volScalarField humiditySource
(
 IOobject
 (
  "humiditySource",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
volScalarField sourceBlanking
(
 IOobject
 (
  "sourceBlanking",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 -1+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
const IOdictionary sourceDict
(
 IOobject
 (
  "sourceDict",    // dictionary name
  runTime.constant(),     // dict is found in "constant"
  mesh,                   // registry for the dict
  IOobject::MUST_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  )
 );

Switch sourceEffects=sourceDict.lookupOrDefault<bool>("sourceEffects", false);
if(sourceEffects)
  {
    List <List<scalar> > sourceInfo(sourceDict.lookup("sourceinfotable"));
    scalarField sourceNumber(sourceInfo.size(),0.0);
    scalarField sourceType(sourceInfo.size(),0.0);
    // Source type 101-199 are velocity source
    // Source type 201-299 are temperature source
    // Source type 301-399 are humidity source
    scalarField sourceTerm1(sourceInfo.size(),0.0);
    scalarField sourceTerm2(sourceInfo.size(),0.0);
    scalarField sourceTerm3(sourceInfo.size(),0.0);
    scalarField sourceTerm4(sourceInfo.size(),0.0);
    scalarField sourceTerm5(sourceInfo.size(),0.0);
    scalarField sourceTerm6(sourceInfo.size(),0.0);
    forAll(sourceInfo,i)
      {
	sourceNumber[i]=sourceInfo[i][0];
	sourceType[i]=sourceInfo[i][1];
	sourceTerm1[i]=sourceInfo[i][2];
	sourceTerm2[i]=sourceInfo[i][3];
	sourceTerm3[i]=sourceInfo[i][4];
	sourceTerm4[i]=sourceInfo[i][5];
	sourceTerm5[i]=sourceInfo[i][6];
	sourceTerm6[i]=sourceInfo[i][7];
      }
    forAll(sourceBlanking,cellI)
      {
	if(sourceBlanking[cellI]>0)
	  {
	    // Term 201 - Heat Release
	    // Term 202 - Exhaust Temperature
	    if(sourceType[sourceBlanking[cellI]-1]==202)
	      temperatureSource[cellI]=sourceTerm1[sourceBlanking[cellI]-1];
	  }
      }
  }

