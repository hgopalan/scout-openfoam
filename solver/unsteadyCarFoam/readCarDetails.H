const IOdictionary carDict
(
 IOobject
 (
  "carDict",    // dictionary name                                                                          
  runTime.constant(),     // dict is found in "constant"                                                                
  mesh,                   // registry for the dict                                                                       
  IOobject::MUST_READ,    // must exist, otherwise failure                                                               
  IOobject::NO_WRITE      // dict is only read by the solver                                                             
  )
 );

List <List<scalar> > carInfo(carDict.lookup("carinfotable"));
scalarField carType(carInfo.size(),0.0);
scalarField carX(carInfo.size(),0.0);
scalarField carY(carInfo.size(),0.0);
scalarField carZ(carInfo.size(),0.0);
scalarField carLx(carInfo.size(),0.0);
scalarField carLy(carInfo.size(),0.0);
scalarField carLz(carInfo.size(),0.0);
scalarField carSpeed(carInfo.size(),0.0);
scalarField carRelease(carInfo.size(),0.0);
scalarField carDirection(carInfo.size(),0.0);
scalarField carTermination(carInfo.size(),0.0);
forAll(carInfo,i)
{
  carX[i]=carInfo[i][0];
  carY[i]=carInfo[i][1];
  carZ[i]=carInfo[i][2];
  carLx[i]=carInfo[i][3];
  carLy[i]=carInfo[i][4];
  carLz[i]=carInfo[i][5];
  carSpeed[i]=carInfo[i][6];
  carRelease[i]=carInfo[i][7];
  carDirection[i]=carInfo[i][8];
  carTermination[i]=carInfo[i][9];
}


