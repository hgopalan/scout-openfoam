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
		      if(cellLabel>0 && mesh.time().value()<=carTermination[i])
			{
			  U[cellLabel][0]=carSpeed[i];
			  U[cellLabel][1]=0;
			  U[cellLabel][2]=0;
			}
		      ymot=ymot+0.25;
		    }
		  xmot=xmot+carDirection[i]*0.25;
		}
	      zloc=zloc+1;
	    }
	  U.correctBoundaryConditions();
	}
    }



