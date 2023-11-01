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
  IOobject::NO_WRITE
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
  IOobject::NO_READ,
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
  IOobject::NO_WRITE
  ),
 1e-10+0*T,
 zeroGradientFvPatchField<scalar>::typeName
 );
volVectorField treeDragForce
(
 IOobject
 (
  "treeDragForce",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 0*U,
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
  IOobject::NO_WRITE
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
  IOobject::NO_WRITE
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
  IOobject::NO_WRITE
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
  IOobject::NO_WRITE
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

const IOdictionary treeBlankingDict
(
 IOobject
 (
  "treeBlankingDict",    // dictionary name
  runTime.constant(),     // dict is found in "constant"
  mesh,                   // registry for the dict
  IOobject::MUST_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  )
 );

Switch treeEffects=treeDict.lookupOrDefault<bool>("treeEffects", false);
scalar treeBlankX1=treeBlankingDict.lookupOrDefault<scalar>("MinX",1e15);
scalar treeBlankY1=treeBlankingDict.lookupOrDefault<scalar>("MinY",1e15);
scalar treeBlankX2=treeBlankingDict.lookupOrDefault<scalar>("MaxX",-1e15);
scalar treeBlankY2=treeBlankingDict.lookupOrDefault<scalar>("MaxY",-1e15);
scalar counter=0;
scalar blankcounter=0;
Info<<"Blanking:"<<treeBlankX1<<" "<<treeBlankY1<<
  treeBlankX2<<" "<<treeBlankY2<<endl;
if(treeEffects)
  {
    // Mark Tree Location 
    solverType solverTypeModel_(solverTypeNames_.get("solverType", ABLDict));
#include "readTreeDetails.H"
#include "treeLoc.H"
    forAll(treeLAD,cellI)
      {
	vector cc=mesh.C()[cellI];
	scalar cx=cc[0];
	scalar cy=cc[1];
	scalar cz=cc[2];
	// Check if blanking blanking cell
	if(cx>treeBlankX1 && cx<treeBlankX2 && cy>treeBlankY1 && cy<treeBlankY2 && treeblanking[cellI]>0)
	  {
	    blankcounter=blankcounter+1;
	    treeblanking[cellI]=-1;
	  }
	if(treeblanking[cellI]>0)
	  {
	    counter=counter+1;
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
	      }
	    treeDragCoeff[cellI]=treeCd[index]*treeLAD[cellI];	      
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
reduce(counter,sumOp<scalar>());
reduce(blankcounter,sumOp<scalar>());
Info<<"Trees marked in "<<counter<<" cells; give or take a few misses"<<endl;
Info<<"Trees blanked in "<<blankcounter<<" cells; give or take a few misses"<<endl;
treeblanking.correctBoundaryConditions();
treeblanking.write();
treeblanking.correctBoundaryConditions();
treeDragCoeff.correctBoundaryConditions();
treeLAD.correctBoundaryConditions();
//treeblanking.write();
treeLAD.write();
treeDragCoeff.write();
treeA.correctBoundaryConditions();
treeE.correctBoundaryConditions();
Info<< "Tree Marking Time = " << runTime.elapsedCpuTime() << " s"
<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
<< nl << endl;

