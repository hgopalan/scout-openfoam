/*const IOdictionary treeLocDict
  (
  IOobject
  (
  "treeLocDict",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE      
  )
  );
  List <List<scalar> > treeLocInfo(treeLocDict.lookup("treeLocTable"));*/
scalar treeCounter=1;
forAll(treeInfo,i)
{
  scalar cellBlank=treeCounter;//treeLocInfo[i][0];
  scalar x=treeCenterX[i];
  scalar y=treeCenterY[i];
  scalar zLower=treeCenterZ[i];
  scalar zUpper=treeCenterZ[i]+treeHeight[i];
  scalar diameter=treeGirth[i];
  scalar z=zLower;
  //  Info<<x<<"  "<<y<<"  "<<zLower<<"  "<<zUpper<<"  "<<diameter<<endl;
  while(z<=zUpper)
    {
      scalar xstart=x-0.5*diameter;
      while(xstart<=(x+0.5*diameter))
	{
	  scalar ystart=y-0.5*diameter;
	  while(ystart<=(y+0.5*diameter))
	    {		  
	      point treeLoc(xstart,ystart,z);
	      label cellLabel=mesh.findCell(treeLoc);
	      if(cellLabel>0)
		{
		  treeblanking[cellLabel]=cellBlank;
		}
	      ystart=ystart+1;
	    }
	  xstart=xstart+1;
	}
      z=z+1;
    }
  treeCounter=treeCounter+1;
}
treeblanking.correctBoundaryConditions();
treeblanking.write();

