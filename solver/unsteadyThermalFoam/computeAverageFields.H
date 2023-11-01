forAll(uavg,levelI)
{
  uavg[levelI]=0.0;
  for(int i=0;i<verticallevels[levelI];++i)
    {
      label celllabel=verticalCellLevels[levelI][i];
      scalar uvalue=U[celllabel][0];
      uavg[levelI]=uavg[levelI]+uvalue;
    }
  uavg[levelI]=uavg[levelI]/(verticallevels[levelI]+SMALL);
  reduce(uavg[levelI],sumOp<scalar>());
  label nproc = Pstream::myProcNo();
  reduce(nproc,maxOp<label>());
  uavg[levelI]=uavg[levelI]/(nproc+1);
}
forAll(vavg,levelI)
{
  vavg[levelI]=0.0;
  for(int i=0;i<verticallevels[levelI];++i)
    {
      label celllabel=verticalCellLevels[levelI][i];
      scalar uvalue=U[celllabel][1];
      vavg[levelI]=vavg[levelI]+uvalue;
    }
  vavg[levelI]=vavg[levelI]/(verticallevels[levelI]+SMALL);
}
/*forAll(Tavg,levelI)
{
  Tavg[levelI]=0.0;
  for(int i=0;i<verticallevels[levelI];++i)
    {
      label celllabel=verticalCellLevels[levelI][i];
      scalar uvalue=T[celllabel];
      Tavg[levelI]=Tavg[levelI]+uvalue;
    }
  Tavg[levelI]=Tavg[levelI]/(verticallevels[levelI]+SMALL);
  }*/

