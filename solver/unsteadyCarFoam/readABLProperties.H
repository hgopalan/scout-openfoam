singlePhaseTransportModel laminarTransport(U, phi);
dimensionedScalar Pr("Pr", dimless, laminarTransport);
dimensionedScalar Prt("Prt", dimless, laminarTransport);
dimensionedScalar rhoRef("rhoRef",dimMass/(dimLength*dimLength*dimLength), laminarTransport);
dimensionedScalar cpRef("cpRef",dimless, laminarTransport);
//dimensionedScalar debugMode("debugMode",dimless,laminarTransport);

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
scalarIOList RHprofile_
(
 IOobject
 (
  "RHProfile",
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
scalarIOList UUpperProfile_
(
 IOobject
 (
  "UUpperProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure                                  
  IOobject::NO_WRITE      // dict is only read by the solver                                
  ),
 dummyscalar
 );
scalarIOList VUpperProfile_
(
 IOobject
 (
  "VUpperProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure                                  
  IOobject::NO_WRITE      // dict is only read by the solver                                
  ),
 dummyscalar
 );
scalarIOList TUpperProfile_
(
 IOobject
 (
  "TUpperProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  ),
 dummyscalar
 );
scalarIOList qVUpperProfile_
(
 IOobject
 (
  "qVUpperProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure                                            
  IOobject::NO_WRITE      // dict is only read by the solver                                         
  ),
 dummyscalar
 );

//Info<<"Profile Size:"<<profilesize<<endl;
List <List<scalar> > profileTable(ABLDict.lookup("profileTable"));
scalarField zProfile(profileTable.size(),0.0);
forAll(zProfile,i)
{
  timeprofile_[i] = profileTable[i][0];
  uprofile_[i] = profileTable[i][1];
  vprofile_[i] = profileTable[i][2];
  Tprofile_[i] = profileTable[i][3];
  RHprofile_[i]= profileTable[i][4];
  MOprofile_[i] = profileTable[i][5];      
  PBLHprofile_[i]=profileTable[i][6];
  MBLProfile_[i]=profileTable[i][7];
  UUpperProfile_[i]=profileTable[i][8];
  VUpperProfile_[i]=profileTable[i][9];
  TUpperProfile_[i]=profileTable[i][10];
  qVUpperProfile_[i]=profileTable[i][11];
}

scalarField timevalues(timeprofile_.size(),1e20);
scalarField uvalues(timeprofile_.size(),0.0);
scalarField vvalues(timeprofile_.size(),0.0);
scalarField rhvalues(timeprofile_.size(),0.0);
scalarField movalues(timeprofile_.size(),0.0);
scalarField pblhvalues(timeprofile_.size(),0.0);
scalarField mblvalues(timeprofile_.size(),0.0);
scalarField TRefvalues(timeprofile_.size(),0.0);
scalarField UUppervalues(timeprofile_.size(),0.0);
scalarField VUppervalues(timeprofile_.size(),0.0);
scalarField TUppervalues(timeprofile_.size(),0.0);
scalarField qVUppervalues(timeprofile_.size(),0.0);
forAll(timeprofile_,i)
{
  timevalues[i] = timeprofile_[i];
  uvalues[i] = uprofile_[i];
  vvalues[i] = vprofile_[i];
  rhvalues[i]= RHprofile_[i];
  movalues[i] = MOprofile_[i];
  pblhvalues[i]=PBLHprofile_[i];
  mblvalues[i]=MBLProfile_[i];
  TRefvalues[i]=Tprofile_[i];
  UUppervalues[i]=UUpperProfile_[i];
  VUppervalues[i]=VUpperProfile_[i];
  TUppervalues[i]=TUpperProfile_[i];
  qVUppervalues[i]=qVUpperProfile_[i];
  //  Info<<"Profile initialized:"<<timevalues[i]<<"  "<<TRefvalues[i]<<endl;
}
//scalar TRef=interpolateXY(mesh.time().value(),timevalues,TRefvalues); 
//scalar MOL=interpolateXY(mesh.time().value(),timevalues,movalues); 
//scalar PBLH=interpolateXY(mesh.time().value(),timevalues,pblhvalues); 
//scalar MBL=interpolateXY(mesh.time().value(),timevalues,mblvalues);
//scalar RH=interpolateXY(mesh.time().value(),timevalues,rhvalues);
//Info<<"Original TRef:"<<TRef<<endl;

