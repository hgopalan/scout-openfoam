forAll(verticallevels,level)
{
  uavg[level]=0.0;
  vavg[level]=0.0;
  verticallevels[level]=0;
}

forAll(U,cellI)
{
  uavg[verticallevelMarker[cellI]]=uavg[verticallevelMarker[cellI]]+U[cellI][0];
  vavg[verticallevelMarker[cellI]]=vavg[verticallevelMarker[cellI]]+U[cellI][1];
  verticallevels[verticallevelMarker[cellI]]=verticallevels[verticallevelMarker[cellI]]+1;
}

forAll(verticallevels,level)
{
  reduce(uavg[level],sumOp<scalar>());
  reduce(vavg[level],sumOp<scalar>());
  reduce(verticallevels[level],sumOp<scalar>());
  uavg[level]=uavg[level]/(verticallevels[level]+1e-5);
  vavg[level]=vavg[level]/(verticallevels[level]+1e-5);
}

