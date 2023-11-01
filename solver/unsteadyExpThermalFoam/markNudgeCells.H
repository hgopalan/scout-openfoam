// Test
//const fvPatchList patches = mesh.boundary();
{
  scalar xWest (ABLDict.lookupOrDefault<scalar>("xWest",-1e30));
  scalar xEast (ABLDict.lookupOrDefault<scalar>("xEast",1e30));
  scalar ySouth (ABLDict.lookupOrDefault<scalar>("ySouth",-1e30));
  scalar yNorth (ABLDict.lookupOrDefault<scalar>("yNorth",1e30));
  scalar zLower (ABLDict.lookupOrDefault<scalar>("zLower",-1e30));
  scalar zUpper (ABLDict.lookupOrDefault<scalar>("zUpper",1e30));
  scalar xCoeff (ABLDict.lookupOrDefault<scalar>("xCoeff",0.0));
  scalar yCoeff (ABLDict.lookupOrDefault<scalar>("yCoeff",0.0));
  scalar zCoeff (ABLDict.lookupOrDefault<scalar>("zCoeff",0.0));
  Info<<"Buffer:"<<endl;
  Info<<xWest<<"  "<<xEast<<"   "<<xCoeff<<endl;
  scalar totalcounter=0;
  scalar counter=0;
  volScalarField cx=(mesh.C()).component(vector::X);
  volScalarField cy=(mesh.C()).component(vector::Y);
  volScalarField cz=(mesh.C()).component(vector::Y);
  scalar domainminx=(min(cx)).value();
  scalar domainmaxx=(max(cx)).value();
  reduce(domainminx,minOp<scalar>());
  reduce(domainmaxx,maxOp<scalar>());
  scalar domainminy=(min(cy)).value();
  scalar domainmaxy=(max(cy)).value();
  reduce(domainminy,minOp<scalar>());
  reduce(domainmaxy,maxOp<scalar>());
  scalar domainminz=(min(cz)).value();
  scalar domainmaxz=(max(cz)).value();
  reduce(domainminz,minOp<scalar>());
  reduce(domainmaxz,maxOp<scalar>());
  scalar coeff=min(xCoeff,yCoeff);
  coeff=min(coeff,zCoeff);
  coeff=min(coeff/3.0,5e-4);
  //scalar C1=0.8;
  //  scalar C2=4;
  /*  scalar meshX=domainmaxx-domainminx;
      scalar meshY=domainmaxy-domainminy;
      scalar xLeft=domainminx+0.2*meshX;
      scalar xRight=domainmaxx-0.2*meshX;
      scalar yLeft=domainminy+0.2*meshY;
      scalar yRight=domainmaxy-0.2*meshY; */
  forAll(blankcell,cellI)
    {
      vector cellCenter=mesh.C()[cellI];
      blankcell[cellI]=0*blankcell[cellI];
      if(cellCenter[2]>=zLower)
	{
	  scalar zBottom=zLower;
	  scalar zTop=zUpper;
	  //	  scalar dampCoeff=(cellCenter[2]-zBottom)/(zTop-zBottom);
	  /*	  if(dampCoeff<0.5)
	    blankcell[cellI]=0.5*coeff*(1-Foam::cos(dampCoeff*3.14159));
	  else
	  blankcell[cellI]=0.5*coeff*(1+(dampCoeff-0.5)*3.14159);      */
	  blankcell[cellI]=coeff*Foam::sqr(Foam::sin(0.5*3.14159*(1-(zTop-cellCenter[2])/(zTop-zBottom))));
	  counter=counter+1;
	}

      if(cellCenter[0]>xEast)
	{
	  scalar zBottom=xEast;
          scalar zTop=domainmaxx;
          scalar dampCoeff=(cellCenter[0]-zBottom)/(zTop-zBottom);
          if(dampCoeff<0.5)
            blankcell[cellI]+=0.5*coeff*(1-Foam::cos(dampCoeff*3.14159));
          else
            blankcell[cellI]+=0.5*coeff*(1+(dampCoeff-0.5)*3.14159);
	  counter=counter+1;
	}
      else if(cellCenter[0]<xWest)
        {
          scalar zBottom=xWest;
          scalar zTop=domainminx;
          scalar dampCoeff=(cellCenter[0]-zBottom)/(zTop-zBottom);
          if(dampCoeff<0.5)
            blankcell[cellI]+=0.5*coeff*(1-Foam::cos(dampCoeff*3.14159));
          else
            blankcell[cellI]+=0.5*coeff*(1+(dampCoeff-0.5)*3.14159);
	  counter=counter+1;
        }     
      if(cellCenter[1]>yNorth)
	{
	  scalar zBottom=yNorth;
          scalar zTop=domainmaxy;
          scalar dampCoeff=(cellCenter[1]-zBottom)/(zTop-zBottom);
          if(dampCoeff<0.5)
            blankcell[cellI]+=0.5*coeff*(1-Foam::cos(dampCoeff*3.14159));
          else
            blankcell[cellI]+=0.5*coeff*(1+(dampCoeff-0.5)*3.14159);
	  counter=counter+1;
	}
      else if(cellCenter[1]<ySouth)
        {
          scalar zBottom=ySouth;
          scalar zTop=domainminy;
          scalar dampCoeff=(cellCenter[1]-zBottom)/(zTop-zBottom);
          if(dampCoeff<0.5)
            blankcell[cellI]+=0.5*coeff*(1-Foam::cos(dampCoeff*3.14159));
          else
            blankcell[cellI]+=0.5*coeff*(1+(dampCoeff-0.5)*3.14159);
	  counter=counter+1;
	  }     
      totalcounter=totalcounter+1;
    }
  reduce(counter, sumOp<scalar>());
  reduce(totalcounter, sumOp<scalar>());
  Info<<"Forced wind speed values:"<<counter<<" Percentage:"<<counter/totalcounter*100.0<<" %"<<endl;
  blankcell.write();
}
