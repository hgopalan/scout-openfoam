forAll(carInfo,i)
{
  scalar xstart=carX[i];
  scalar ystart=carY[i];
  scalar zstart=carZ[i];
  scalar xend=carX[i]+carDirection[i]*carLx[i];
  scalar yend=carY[i]+carLy[i];
  scalar zend=carZ[i]+carLz[i];
  scalar zloc=zstart;
  scalar release=carRelease[i];
  if(mesh.time().value()>carRelease[i])
    {
      while(zloc<=zend)
	{
	  scalar xmot=xstart;
	  while((xmot<=xend && carDirection[i]==1) || (xmot>=xend && carDirection[i]==-1))
	    {
	      scalar ymot=ystart;
	      while(ymot<=yend)
		{		  
		  point carLoc(xmot,ymot,zloc);
		  label cellLabel=mesh.findCell(carLoc);
		  if(cellLabel>0)
		    {
		      PTraffic[cellLabel]=Foam::sqrt(0.5*0.9*Foam::pow(carSpeed[i],3)*carLy[i]);
		    }
		  ymot=ymot+0.25;
		}
	      xmot=xmot+carDirection[i]*0.25;
	    }
	  zloc=zloc+1;
	}
    }
}



