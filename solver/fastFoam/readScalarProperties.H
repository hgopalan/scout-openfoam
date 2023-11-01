const IOdictionary scalarDict
(
 IOobject
 (
  "scalarDict",    // dictionary name
  runTime.constant(),     // dict is found in "constant"
  mesh,                   // registry for the dict
  IOobject::MUST_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  )
 );
// Testing code for unsteady inflow 
Switch runScalar(scalarDict.lookupOrDefault<bool>("runScalar",false));
scalar scalarsize_ (scalarDict.lookupOrDefault<scalar>("scalarSize",1));
scalarField dummyscalar1(scalarsize_,0.0);
scalarIOList gasDensity_
(
 IOobject
 (
  "gasDensity",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummyscalar1
 ); 
scalarIOList initialRelTime_
(
 IOobject
 (
  "initialRelTime",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummyscalar1
 ); 
scalarIOList finalRelTime_
(
 IOobject
 (
  "finalRelTime",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummyscalar1
 ); 


List <List<scalar> > scalarTable(scalarDict.lookup("profileTable"));
scalarField zProfile1(scalarTable.size(),0.0);
forAll(zProfile1,i)
{
  gasDensity_[i] = scalarTable[i][0];
  initialRelTime_[i] = scalarTable[i][1];
  finalRelTime_[i] = scalarTable[i][2];
}
Info<<"Done reading scalar properties for:"<<scalarsize_<<" scalars"<<endl;


