scalarField zLoc(1024,1);
scalarField vertT(1024,TRef);
scalarField vertqV(1024,1e-4);
scalarField vertNut(1024,0.1);
scalarField vertU(1024,1.0);
scalarField vertV(1024,0.0);
scalar lapseRate=0.004;
scalar zLower=0;
scalar zUpper=0;
scalar upperMarked=0;
scalar lowerMarked=0;
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
forAll(zLoc,i)
{
  zLoc[i]=i*2;
}
Info<<"Creating vertical profiles for Temperature and mixing ratio"<<endl;
