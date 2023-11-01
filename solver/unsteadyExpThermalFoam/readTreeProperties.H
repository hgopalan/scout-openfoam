// Initialize tree fields 
// Tree-lead Temperature 
Random perturbation(1234567);
volScalarField leafTemperature
(
 IOobject
 (
  "leafTemperature",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
// Boolean variable for marking trees 
volScalarField treeblanking
(
 IOobject
 (
  "treeBlanking",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 -1+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
// Leaf area density variation of the tree 
volScalarField treeLAD
(
 IOobject
 (
  "treeLAD",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 1e-10+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
// Drag Force due to the tree 
volScalarField treeDragCoeff
(
 IOobject
 (
  "treeDragCoeff",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 1e-10+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
// Shading of the tree 
volScalarField treeShading
(
 IOobject
 (
  "treeShading",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 1+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
// Absoptivity 
volScalarField treeA
(
 IOobject
 (
  "treeA",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 1e-10+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
volScalarField treeE
(
 IOobject
 (
  "treeE",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 1e-10+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
volScalarField treeTau
(
 IOobject
 (
  "treeTau",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 1e-10+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );

const IOdictionary treeDict
(
 IOobject
 (
  "treeDict",    // dictionary name
  runTime.constant(),     // dict is found in "constant"
  mesh,                   // registry for the dict
  IOobject::MUST_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  )
 );

Switch treeEffects=treeDict.lookupOrDefault<bool>("treeEffects", false);
if(treeEffects)
  {
    // Dummy 
    // Tree is characterized by the following parameters 
    // Type: (1) Shrub or Forest ; (2) Individual tree ; (3) Hanging tree/shrub 
    // Height (in meters)
    // LAI 
    // Coefficient of drag Cd 
    // Zm - Location of highest LAI for model 2; Start of LAI for model 3 
    // Lm - Calculated for individual tree model 
    // Absorptivity 
    // Emissivity 
    List <List<scalar> > treeInfo(treeDict.lookup("treeinfotable"));
    scalarField treeType(treeInfo.size(),0.0);
    scalarField treeHeight(treeInfo.size(),0.0);
    scalarField treeLAI(treeInfo.size(),0.0);
    scalarField treeCd(treeInfo.size(),0.0);
    scalarField treeZm(treeInfo.size(),0.0);
    scalarField treeLm(treeInfo.size(),0.0);
    scalarField treeAbsorptivity(treeInfo.size(),0.0);
    scalarField treeEmissivity(treeInfo.size(),0.0);
    scalarField treeTransmissivity(treeInfo.size(),0.0);
    forAll(treeInfo,i)
      {
	treeType[i]=treeInfo[i][0];
	treeHeight[i]=treeInfo[i][1];
	treeLAI[i]=treeInfo[i][2];
	treeCd[i]=treeInfo[i][3];
	treeZm[i]=treeInfo[i][4]*treeHeight[i];
	treeAbsorptivity[i]=treeInfo[i][5];
	treeEmissivity[i]=treeInfo[i][6];
	treeTransmissivity[i]=treeInfo[i][7];
	if(treeType[i]==2)
	  {
	    scalar z=0;
	    scalar dz=treeHeight[i]/100;
	    scalar expFun=0;
	    while(z<=treeHeight[i])
	      {
		scalar ratio=(treeHeight[i]-treeZm[i])/(treeHeight[i]-z+SMALL);
		if(z<treeZm[i])
		  {
		    expFun=expFun+Foam::pow(ratio,6.0)*Foam::exp(6*(1-ratio));
		  }
		else
		  {
		    expFun=expFun+Foam::pow(ratio,0.5)*Foam::exp(0.5*(1-ratio));
		  }
		z=z+dz;
	      }
	    treeLm[i]=treeLAI[i]/(expFun*dz);
	    //	    Info<<"Maximum LM:"<<treeLm[i]<<endl;
	  }
      }
    forAll(treeLAD,cellI)
      {
	vector cc=mesh.C()[cellI];
	scalar cz=cc[2];
	if(treeblanking[cellI]>0)
	  {
	    scalar index=int(treeblanking[cellI])-1;
	    if(treeType[index]==1)
	      {
		scalar eta=min(cz/treeHeight[index],1.0);
		treeLAD[cellI]=treeLAI[index]/treeHeight[index]*Foam::exp(-10*(1-eta));
	      }
	    else if(treeType[index]==2)
	      {
		scalar ratio=(treeHeight[index]-treeZm[index])/(treeHeight[index]-cz+1e-10);
		if(cz<treeZm[index])
		  treeLAD[cellI]=treeLm[index]*Foam::pow(ratio,6.0)*Foam::exp(6*(1-ratio));
		else if(cz<=treeHeight[index])
		  treeLAD[cellI]=treeLm[index]*Foam::pow(ratio,0.5)*Foam::exp(0.5*(1-ratio));
		//		Pout<<cz<<"  "<<treeLAD[cellI]<<"   "<<treeZm[index]<<"   "<<treeLm[index]<<"  "<<treeHeight[index]<<endl;
	      }
	    //	    else if(treeType[index]==3)
	    //	      treeLAD[cellI]=pos(cz-treeZm[index])*treeLAI[index]/(treeHeight[index]-treeZm[index]);
	    treeDragCoeff[cellI]=treeCd[index]*treeLAD[cellI];	      
	    //	    scalar fmu=sqr(1-cz/treeHeight[index]);
	    treeShading[cellI]=treeHeight[index];
	    treeA[cellI]=treeLAD[cellI]*treeAbsorptivity[index]*Foam::pow(mesh.V()[cellI],1.0/3.0);
	    treeE[cellI]=treeLAD[cellI]*treeEmissivity[index]*Foam::pow(mesh.V()[cellI],1.0/3.0);
	    treeTau[cellI]=treeLAD[cellI]*treeTransmissivity[index]*Foam::pow(mesh.V()[cellI],1.0/3.0);
	  }
      }
  }
//treeA*=treeLAD/gMax(treeLAD);
//treeE*=treeLAD/gMax(treeLAD);
//treeTau*=treeLAD/gMax(treeLAD);
treeDragCoeff.correctBoundaryConditions();
treeLAD.correctBoundaryConditions();
treeblanking.correctBoundaryConditions();
treeShading.correctBoundaryConditions();
treeLAD.write();
treeA.correctBoundaryConditions();
treeE.correctBoundaryConditions();
