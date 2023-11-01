    // Mark Building 
    buildingBlanking=0*buildingBlanking+1;
    forAll(buildingBlanking,cellI)
      {
	if(regionBlanking[cellI]>0 && regionBlanking[cellI]<=regionCounter)
	  {
	    buildingBlanking[cellI]=0;
	    treeA[cellI]=0.8;
	    treeE[cellI]=0.9;
	  }
      }
    buildingBlanking.correctBoundaryConditions();
    buildingBlankingSurface=fvc::interpolate(buildingBlanking);
    buildingBlankingSurface=pos(buildingBlankingSurface-0.8)*buildingBlankingSurface; 
	treeA.correctBoundaryConditions();
    treeE.correctBoundaryConditions();
	// Mark Cells 
	List<label> blankCellList;
	List<scalar> internalTemperature;
    forAll(buildingBlanking,cellI)
      {
	if(regionBlanking[cellI]>0 && regionBlanking[cellI]<=regionCounter)
	  {
		blankCellList.append(cellI);
		scalar temperature=298+0.01*mesh.C()[cellI][2]+0.01*Foam::mag(mesh.C()[cellI][0])+0.01*Foam::mag(mesh.C()[cellI][1]);
		internalTemperature.append(temperature);
	  }
      }
