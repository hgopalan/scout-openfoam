PtrList<volVectorField> Urain;
PtrList<surfaceScalarField> phirain;
PtrList<volScalarField> alpharain;
PtrList<volScalarField> scr;
volScalarField blankRain=max(neg(treeblanking),0.5+0*treeblanking);
Urain.setSize(phases.size());
phirain.setSize(phases.size());
alpharain.setSize(phases.size());
scr.setSize(phases.size());

std::stringstream converter; 	
for (label phase_no = 0; phase_no < phases.size(); phase_no++)
  {
    string phase_Uname = "U";
    string phase_phiname = "phi";
    string phase_alphaname = "alpha";
	
    string phase_no_str;
    std::stringstream convert;
    convert << phase_no+1;
    phase_no_str = convert.str();
	
    phase_Uname.append(phase_no_str);
    phase_phiname.append(phase_no_str);
    phase_alphaname.append(phase_no_str);

    Urain.set
      (
       phase_no,
       new volVectorField
       (
	IOobject
	(
	 phase_Uname,
	 "0",
	 //	runTime.timeName(),
	 mesh,
	 IOobject::MUST_READ,
	 IOobject::AUTO_WRITE
	 ),
	mesh
	)
       );
    //	Urain[phase_no] = Utemp;

    phirain.set
      (
       phase_no,
       new surfaceScalarField
       (
	IOobject
	(
	 phase_phiname,
	 "0",
	 //	 runTime.timeName(),
	 mesh,
	 IOobject::NO_READ,
	 IOobject::NO_WRITE
	 ),
	linearInterpolate(Urain[phase_no]) & mesh.Sf()
	)
       );
    //	phirain[phase_no] = phitemp;

    alpharain.set
      (
       phase_no,
       new volScalarField 
       (
	IOobject
	(
	 phase_alphaname,
	 "0",
	 //	 runTime.timeName(),
	 mesh,
	 IOobject::MUST_READ,
	 IOobject::AUTO_WRITE
	 ),
	mesh
	)
       );
    //	alpharain[phase_no] = alphatemp;
	
  }
// Cheat U 
volVectorField UCheat
(
 IOobject
 (
  "UCheat",
  "0",
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 U,
 zeroGradientFvPatchField<vector>::typeName
 );
UCheat=U;
UCheat.correctBoundaryConditions();	
