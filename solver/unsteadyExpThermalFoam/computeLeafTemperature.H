{
  //  const objectRegistry& db = mesh.thisDb();
  //  const volScalarField& totalRadiation=db.lookupObject<volScalarField>("G");
  //  const volScalarField& swRadiation=db.lookupObject<volScalarField>("qcb");
  // Aerodynamic Resistance
  volScalarField qsenLeaf=0*T;
  volScalarField swRad=mag(fvc::grad(swRadiation));
  volScalarField lwRad=mag(fvc::grad(totalRadiation))-swRad;
  forAll(leafTemperature,cellI)
    {
      scalar eps=100;
      scalar counter=0;
      scalar fx=0;
      scalar fxp=0;
      scalar radiation=0.0;
      if(treeblanking[cellI]>0)
	{
	  scalar ra_=130*Foam::sqrt(0.5/(mag(U[cellI])+0.05));
	  while(counter<10)
	    {
	      radiation=swRad[cellI]*(treeA[cellI]-treeTau[cellI])+lwRad[cellI]*treeE[cellI];
	      fx=rhok[cellI]*cpRef.value()/ra_*(leafTemperature[cellI]-T[cellI]);
	      fx=fx-radiation/treeLAD[cellI];
	      fx=fx+treeE[cellI]*5.67e-8*pow(leafTemperature[cellI],4)-treeE[cellI]*(213+5.5*(T[cellI]-273.15));
	      fxp=rhok[cellI]*cpRef.value()/ra_+4*treeE[cellI]*5.67e-8*pow(leafTemperature[cellI],3);
	      eps=fx/(fxp+SMALL);
	      leafTemperature[cellI]=leafTemperature[cellI]-eps;
	      counter=counter+1;
	    }
	  //	  Pout<<fx/fxp<<"   "<<ra_<<"   "<<treeLAD[cellI]<<"   "<<lwRad[cellI]<<"  "
	  //	      <<swRad[cellI]<<"   "<<radiation<<"   "<<leafTemperature[cellI]<<"   "<<T[cellI]<<endl;
	  qsenLeaf[cellI]=2*rhok[cellI]*cpRef.value()/ra_*(leafTemperature[cellI]-T[cellI]);
	}
    }
  leafTemperature.correctBoundaryConditions();
  sourceT=treeLAD/(rhok*cpRef)*qsenLeaf;
}
