// Dummy 
// Tree is characterized by the following parameters 
// Type: (1) Shrub or Forest ; (2) Individual tree ; (3) Hanging tree/shrub 
// Height (in meters)
// LAI 
// Coefficient of drag Cd 
// Zm - Location of highest LAI for model 2; Start of LAI for model 3 
// Lm - Calculated for individual tree model 
// Absorptivity 
// Emissivity 
List <List<scalar> > treeInfo(treeDict.lookup("treeinfotable"));
scalarField treeType(treeInfo.size(),0.0);
scalarField treeCenterX(treeInfo.size(),0.0);
scalarField treeCenterY(treeInfo.size(),0.0);  
scalarField treeCenterZ(treeInfo.size(),0.0);  
scalarField treeGirth(treeInfo.size(),0.0);    
scalarField treeHeight(treeInfo.size(),0.0);
scalarField treeLAI(treeInfo.size(),0.0);
scalarField treeCd(treeInfo.size(),0.0);
scalarField treeZm(treeInfo.size(),0.0);
scalarField treeLm(treeInfo.size(),0.0);
scalarField treeAbsorptivity(treeInfo.size(),0.0);
scalarField treeEmissivity(treeInfo.size(),0.0);
scalarField treeTransmissivity(treeInfo.size(),0.0);
forAll(treeInfo,i)
{
  treeType[i]=treeInfo[i][0];
  treeCenterX[i]=treeInfo[i][1];
  treeCenterY[i]=treeInfo[i][2];
  treeCenterZ[i]=treeInfo[i][3];
  treeHeight[i]=treeInfo[i][4];
  treeGirth[i]=treeInfo[i][5];	
  treeLAI[i]=treeInfo[i][6];
  treeCd[i]=treeInfo[i][7];
  treeZm[i]=treeInfo[i][8]*treeHeight[i];
  treeAbsorptivity[i]=treeInfo[i][9];
  treeEmissivity[i]=treeInfo[i][10];
  treeTransmissivity[i]=treeInfo[i][11];
  // If the tree is small need to switch model and LAI
  if(treeHeight[i]<10.0)
    {
      treeType[i]=1;
      treeLAI[i]=1.0;
      treeCd[i]=0.1;
    }

  switch (solverTypeModel_)
    {
    case steadyWind:
      {
	treeAbsorptivity[i]=0;
	treeEmissivity[i]=0;
	treeTransmissivity[i]=0;
	break;
      }
    case unsteadyWind:
      {
	treeAbsorptivity[i]=0;
	treeEmissivity[i]=0;
	treeTransmissivity[i]=0;
	break;
      } 
    case steadyThermal:
      break;
    case unsteadyThermal:
      break;
    case steadyHumid:
      break;
    case unsteadyHumid:
      break;
    }
  if(treeType[i]==2)
    {
      scalar z=0;
      scalar dz=treeHeight[i]/100;
      scalar expFun=0;
      while(z<=treeHeight[i])
	{
	  scalar ratio=(treeHeight[i]-treeZm[i])/(treeHeight[i]-z+SMALL);
	  if(z<treeZm[i])
	    {
	      expFun=expFun+Foam::pow(ratio,6.0)*Foam::exp(6*(1-ratio));
	    }
	  else
	    {
	      expFun=expFun+Foam::pow(ratio,0.5)*Foam::exp(0.5*(1-ratio));
	    }
	  z=z+dz;
	}
      treeLm[i]=treeLAI[i]/(expFun*dz);
      //	    Info<<"Maximum LM:"<<treeLm[i]<<endl;
    }
}

