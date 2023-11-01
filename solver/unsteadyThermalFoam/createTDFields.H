PtrList<volScalarField> Ctrain;
Ctrain.setSize(phases.size());	
for (label phase_no = 0; phase_no < phases.size(); phase_no++)
  {
    string phase_Ctname = "Ct";
    string phase_no_str;         
    std::stringstream convert;
    convert << phase_no+1;
    phase_no_str = convert.str();
    phase_Ctname.append(phase_no_str);

    Ctrain.set
      (
       phase_no,
       new volScalarField
       (
	IOobject
	(
	 phase_Ctname,
	 "0",
	 mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::NO_WRITE
	 ),
	mesh,
	dimensionedScalar(phase_Ctname,dimensionSet(0,0,0,0,0,0,0),0)
	)
       );
    //    Ctrain.set(phase_no,Cttemp);
  }

volScalarField nutrain
(
 IOobject
 (
  "nutrain",
  "0",
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("nutrain",dimensionSet(0,2,-1,0,0,0,0),1)
 );

volScalarField tp
(
 IOobject
 (
  "tp",
  "0",
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("tp",dimensionSet(0,0,1,0,0,0,0),0)
 );

volScalarField tfl
(
 IOobject
 (
  "tfl",
  "0",
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("tfl",dimensionSet(0,0,1,0,0,0,0),0)
 );
