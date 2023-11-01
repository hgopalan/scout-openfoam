scalarField zLoc(1300,1);
scalarField vertT(1300,TRef);
scalarField vertqV(1300,1e-4);
scalarField vertNut(1300,0.1);
scalarField vertU(1300,1.0);
scalarField vertV(1300,0.0);
scalarIOList zProfile_
(
 IOobject
 (
  "zProfile",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 zLoc
 ); 
scalarIOList vertUProfile_
(
 IOobject
 (
  "vertUProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure                          
  IOobject::NO_WRITE      // dict is only read by the solver                                         
  ),
 zLoc
 );
scalarIOList vertVProfile_
(
 IOobject
 (
  "vertVProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure                                         
  IOobject::NO_WRITE      // dict is only read by the solver                                         
  ),
 zLoc
 );

scalarIOList vertNutProfile_
(
 IOobject
 (
  "vertNutProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure                                        
  IOobject::NO_WRITE      // dict is only read by the solver                                         
  ),
 zLoc
 );
 scalarIOList vertTProfile_
(
 IOobject
 (
  "vertTProfile",
  runTime.constant(),  
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver  
  ),
 zLoc
 );
scalarIOList vertqVProfile_
(
 IOobject
 (
  "vertqVProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver                                                                  
  ),
 zLoc
 );
//scalar extraDepth=(3000-200)/150;
forAll(zLoc,i)
{
  // if(i<=99)
  zLoc[i]=i*2;
  // else
  //   zLoc[i]=zLoc[i-1]+extraDepth;
}
Info<<"Creating vertical profiles for Temperature and mixing ratio"<<endl;
