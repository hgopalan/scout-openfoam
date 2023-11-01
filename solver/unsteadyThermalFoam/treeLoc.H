scalar treeCounter=0;
//scalar partialTreeCounter=0;
//scalar PointCounter=0;
//scalar totalTrees=treeInfo.size();
//label procID=Pstream::myProcNo();
forAll(treeInfo,i)
{
  scalar x=treeCenterX[i];
  scalar y=treeCenterY[i];
  // Break corners 
  scalar radius=0.65*treeGirth[i];
  scalar zUpper=treeCenterZ[i]+treeHeight[i];
  //PointCounter=0;
  forAll(U,cellI)
    {
      vector cell=mesh.C()[cellI];
      scalar cx=cell[0];
      scalar cy=cell[1];
      scalar cz=cell[2];
      // 1.1 to break corners 
      scalar distance=Foam::sqrt(Foam::pow(cx-x,2)+Foam::pow(cy-y,2));
      if(distance<=radius && cz>=treeCenterZ[i] && cz<=zUpper)
	{
	  treeblanking[cellI]=i;
	  //PointCounter=PointCounter+1;
	}
    }
  scalar maxPoints=1.2*treeGirth[i]*treeGirth[i]*treeHeight[i];
  /*if(PointCounter>0 && PointCounter< (0.5*maxPoints))
    partialTreeCounter=partialTreeCounter+1;
  else if(PointCounter>=(0.5*maxPoints))
  treeCounter=treeCounter+1;*/
  if(i%100==0)
    Info<<"Trees progress:"<<i<<" out of:"<<treeInfo.size()<<endl;
}
//Pout<<"Processor:"<<procID<<" Partial Trees:"<<partialTreeCounter<<
// " Full Trees:"<<treeCounter<<" out of "<<totalTrees<<endl;
reduce(treeCounter,sumOp<scalar>());
Info<<"Trees Done"<<endl;
