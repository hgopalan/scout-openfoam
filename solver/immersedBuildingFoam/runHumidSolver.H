Info<<" Starting Humid Solver"<<endl;
#include "adjustTurbulenceConstants.H"
turb_C3.write();
//#include "runPotentialFlowCorrection.H"
//if(mesh.time().value()==0)
//  radiation->correct();
blankcell.write();
vector vertical(0,0,1.0);
scalar fastCounter=0;
if(hotSpot)
  {
#include "runMassConsistentCorrection.H"
  }
scalar regionCounter=1;
while (runTime.run())
  {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
    runTime++;
    Info<< "Time = " << runTime.timeName() << nl << endl;
    // --- Pressure-velocity PIMPLE corrector loop
    //#include "createNudgeField.H"
    #include "updateICVariables.H"
    #include "updateBuildingLocation.H"
#include "1dInitalize.H" // Create vertical scalar list for Temperature and mixing ratio 
#include "createWindandTurbulence.H" // Create the three-dimensional fields 
#include "computeUstarThetastar.H" // Compute friction coefficients 
#include "1dProfileTemp.H" // One-dimensional vertical Temperature profile
#include "1dProfileqV.H" // One-dimensional mixing ratio profile 
#include "createTemperatureandHumidity.H"    
#include "adjustTurbulenceConstants.H"
    scalar radiationcounter=0;
    
    radiation->correct();
    while (pimple.loop())
      {
	/*List<label> blankGhostCells;
	List<vector> blankGhostWind;
	// Ghost 
	for( label i=0; i < mesh.nInternalFaces(); i++)
	  {
	    if( (buildingBlanking[own[i]]==1 && buildingBlanking[nei[i]]==0)) 
	      {
		U[nei[i]]=-U[own[i]];
		blankGhostCells.append(nei[i]);
		blankGhostWind.append(U[nei[i]]);
	      }
	    if( (buildingBlanking[own[i]]==0 && buildingBlanking[nei[i]]==1)) 
	      {
		U[own[i]]=-U[nei[i]];
		blankGhostCells.append(own[i]);
		blankGhostWind.append(U[own[i]]);
	      }
	  }*/
	tempcounter=0;
	forAll(nudgeU_,cellI)
	  {
	    nudgeU_[cellI][0]=blankcell[cellI]*(Unudge[cellI][0] - U[cellI][0]);
	    nudgeU_[cellI][1]=blankcell[cellI]*(Unudge[cellI][1] - U[cellI][1]);
	    nudgeU_[cellI][2]=blankcell[cellI]*(Unudge[cellI][2] - U[cellI][2]);
	  }
	nudgeNut_=blankcell*(nutnudge - nut);
	// Thermal part comes here
	// Steady we need to cheat
	switch (solverTypeModel_)
	  {
          case steadyWind:
            {
              radiationcounter=0;
              break;
            }
          case unsteadyWind:
            {
              radiationcounter=0;
              break;
            }
          case steadyThermal:
            {
              radiationcounter=0;
              break;
            }
          case unsteadyThermal:
            {
              radiationcounter=0;
              break;
            }
	    
	  case steadyHumid:
	    {
	      radiationcounter=radiationcounter+1;
	      if(radiationcounter == callRadiation)
		{
		  radiation->correct();
		  radiationcounter=0;
		}
	      break;
	    }
	  case unsteadyHumid:
	    {
	      // Unsteady Radiation calculation at final iteration of pimple
	      radiationcounter=radiationcounter+1;
	      if(pimple.finalIter() || radiationcounter==callRadiation)
		{
		  radiation->correct();
		  radiationcounter=0;
		}
	      break;
	    }
	  }
#include "TEqn.H"
#include "qvEqn.H"
#include "computeEvapotranspiration.H"
#include "pmvcalculator.H"
        if(fullSolver)
          {
#include "UEqn.H"
#include"nutEqn.H"
#include "pEqn.H"
	  }      
	if(hotSpot)
	  {
	    U.correctBoundaryConditions();
	    nut.correctBoundaryConditions();
	  }
	fastCounter=fastCounter+1;
      }
    regionCounter=regionCounter+1;      
    runTime.write();
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
	<< nl << endl;
  }
	

