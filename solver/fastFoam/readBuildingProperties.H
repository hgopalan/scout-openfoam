scalar totalBuildings (ABLDict.lookupOrDefault<scalar>("totalBuildings",100));
scalarField dummybuildings(totalBuildings,0.0);
scalarIOList buildingID
(
 IOobject
 (
  "buildingID",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummybuildings
 ); 
 scalarIOList buildingminX
(
 IOobject
 (
  "buildingminX",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummybuildings
 ); 
  scalarIOList buildingmaxX
(
 IOobject
 (
  "buildingmaxX",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummybuildings
 ); 
 scalarIOList buildingminY
(
 IOobject
 (
  "buildingminY",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummybuildings
 ); 
  scalarIOList buildingmaxY
(
 IOobject
 (
  "buildingmaxY",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummybuildings
 );  
 scalarIOList buildingminZ
(
 IOobject
 (
  "buildingminZ",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummybuildings
 ); 
  scalarIOList buildingmaxZ
(
 IOobject
 (
  "buildingmaxZ",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 dummybuildings
 ); 

Info<<"Number of Buildings:"<<totalBuildings<<endl;
List <List<scalar> > buildingprofileTable(ABLDict.lookup("buildingInformation"));
scalarField zzProfile(buildingprofileTable.size(),0.0);
forAll(zzProfile,i)
{
 buildingID[i] = buildingprofileTable[i][0];
 buildingminX[i]=buildingprofileTable[i][1];
 buildingmaxX[i]=buildingprofileTable[i][2];
 buildingminY[i]=buildingprofileTable[i][3];
 buildingmaxY[i]=buildingprofileTable[i][4];
 buildingminZ[i]=buildingprofileTable[i][5];
 buildingmaxZ[i]=buildingprofileTable[i][6];
}
