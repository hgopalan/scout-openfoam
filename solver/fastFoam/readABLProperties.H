singlePhaseTransportModel laminarTransport(U, phi);
dimensionedScalar Pr("Pr", dimless, laminarTransport);
dimensionedScalar Prt("Prt", dimless, laminarTransport);
dimensionedScalar rhoRef("rhoRef",dimMass/(dimLength*dimLength*dimLength), laminarTransport);
dimensionedScalar cpRef("cpRef",dimless, laminarTransport);

const IOdictionary ABLDict
(
 IOobject
 (
  "ABLDict",    // dictionary name
  runTime.constant(),     // dict is found in "constant"
  mesh,                   // registry for the dict
  IOobject::MUST_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  )
 );
//word radMethod(ABLDict.lookup("radiationMethod"));
// Testing code for unsteady inflow 
scalar profilesize (ABLDict.lookupOrDefault<scalar>("profileSize",100));
scalar TZRef_(ABLDict.lookupOrDefault<scalar>("TZRef",2.0));
scalar ZREF_(ABLDict.lookupOrDefault<scalar>("ZRef",10.0));
scalar Z0_(ABLDict.lookupOrDefault<scalar>("z0",1e-4));
Switch TEqn(ABLDict.lookupOrDefault<Switch>("TEqn",true));
Switch Steady(ABLDict.lookupOrDefault<Switch>("Steady",true));
Switch Nudge(ABLDict.lookupOrDefault<Switch>("Nudge",false));
//scalar timeScale_(ABLDict.lookupOrDefault<scalar>("timeScale",3600));
//scalar callRadiation(ABLDict.lookupOrDefault<scalar>("timeScale",50));
Switch multiSteady(ABLDict.lookupOrDefault<Switch>("multiSteady",false));
Switch writeSummary(ABLDict.lookupOrDefault<Switch>("writeSummary",true));
scalarField dummyscalar(profilesize,0.0);
scalarIOList timeprofile_
(
 IOobject
 (
  "timeProfile",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummyscalar
 ); 
scalarIOList uprofile_
(
 IOobject
 (
  "uProfile",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummyscalar
 ); 
scalarIOList vprofile_
(
 IOobject
 (
  "vProfile",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummyscalar
 ); 
scalarIOList Tprofile_
(
 IOobject
 (
  "TProfile",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummyscalar
 );   
scalarIOList MOprofile_
(
 IOobject
 (
  "MOProfile",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummyscalar
 );   
scalarIOList PBLHprofile_
(
 IOobject
 (
  "PBLHProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  ),
 dummyscalar
 );
scalarIOList MBLProfile_
(
 IOobject
 (
  "MBLProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure                                             
  IOobject::NO_WRITE      // dict is only read by the solver                                          
  ),
 dummyscalar
 );
scalarIOList switchProfile_
(
 IOobject
 (
  "switchProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  ),
 dummyscalar
 );

Info<<"Profile Size:"<<profilesize<<endl;
List <List<scalar> > profileTable(ABLDict.lookup("profileTable"));
scalarField zProfile(profileTable.size(),0.0);
forAll(zProfile,i)
{
  timeprofile_[i] = profileTable[i][0];
  uprofile_[i] = profileTable[i][1];
  vprofile_[i] = profileTable[i][2];
  Tprofile_[i] = profileTable[i][3];
  MOprofile_[i] = profileTable[i][4];      
  PBLHprofile_[i]=profileTable[i][5];
  MBLProfile_[i]=profileTable[i][6];
  if(Steady)
    switchProfile_[i]=0.0;
  else
    switchProfile_[i]=1.0;
}

scalarField timevalues(timeprofile_.size(),1e20);
scalarField uvalues(timeprofile_.size(),0.0);
scalarField vvalues(timeprofile_.size(),0.0);
scalarField movalues(timeprofile_.size(),0.0);
scalarField pblhvalues(timeprofile_.size(),0.0);
scalarField mblvalues(timeprofile_.size(),0.0);
scalarField TRefvalues(timeprofile_.size(),0.0);
scalarField switchProfilevalues(timeprofile_.size(),0.0);
forAll(timeprofile_,i)
{
  timevalues[i] = timeprofile_[i];
  uvalues[i] = uprofile_[i];
  vvalues[i] = vprofile_[i];
  movalues[i] = MOprofile_[i];
  pblhvalues[i]=PBLHprofile_[i];
  mblvalues[i]=MBLProfile_[i];
  TRefvalues[i]=Tprofile_[i];
  switchProfilevalues[i]=switchProfile_[i];
  Info<<"Profile initialized:"<<timevalues[i]<<"  "<<TRefvalues[i]<<endl;
}
scalar TRef=interpolateXY(mesh.time().value(),timevalues,TRefvalues); 
scalar MOL=interpolateXY(mesh.time().value(),timevalues,movalues); 
scalar PBLH=interpolateXY(mesh.time().value(),timevalues,pblhvalues); 
scalar MBL=interpolateXY(mesh.time().value(),timevalues,mblvalues); 
//Info<<"Original TRef:"<<TRef<<endl;
// Initialize Buildings
#include "readBuildingProperties.H"
