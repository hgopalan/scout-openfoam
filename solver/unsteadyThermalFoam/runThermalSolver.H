Info<<" Starting Humid Solver"<<endl;
#include "adjustTurbulenceConstants.H"
Switch fullSolver=ABLDict.lookupOrDefault<bool>("fullSolver", true);
Switch hotSpot=ABLDict.lookupOrDefault<bool>("hotSpot", false);
//turb_C3.write();
if(mesh.time().value()==0)
  radiation->correct();
//blankcell.write();
vector vertical(0,0,1.0);
scalar fastCounter=0;
Switch firstTime=true;
if(hotSpot)
  {
#include "runMassConsistentCorrection.H"
  }
while (runTime.run())
  {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
    runTime++;
    Info<< "Time = " << runTime.timeName() << nl << endl;
    // --- Pressure-velocity PIMPLE corrector loop
    //#include "createNudgeField.H"
    locU1=interpolateXY(U.time().value(),timevalues,uvalues);
    locV1=interpolateXY(U.time().value(),timevalues,vvalues);
    M=Foam::sqrt(sqr(locU1)+sqr(locV1));
    vector locFlowDir(locU1/M,locV1/M,0);
    RH=interpolateXY(mesh.time().value(),timevalues,rhvalues);
    TRef=interpolateXY(U.time().value(),timevalues,TRefvalues);
    UUpper=interpolateXY(mesh.time().value(),timevalues,UUppervalues);
    VUpper=interpolateXY(mesh.time().value(),timevalues,VUppervalues);
    TUpper=min(interpolateXY(mesh.time().value(),timevalues,TUppervalues),306);    
    qVUpper=interpolateXY(mesh.time().value(),timevalues,qVUppervalues);
    TRef=interpolateXY(U.time().value(),timevalues,TRefvalues);
    MOL=interpolateXY(U.time().value(),timevalues,movalues);
    if(Foam::mag(MOL)>25)
      MOL=MOL;
    else if(Foam::neg(MOL) && MOL> -25)
      MOL=-25;
    else if(Foam::pos(MOL) && Foam::mag(MOL)<25)
      MOL=25;
    Info<<"Monin-Obukhov Length:"<<MOL<<endl;
    PBLH=interpolateXY(U.time().value(),timevalues,pblhvalues);
    MBL=interpolateXY(U.time().value(),timevalues,mblvalues);
    //#include "1dInitalize.H" // Create vertical scalar list for Temperature and mixing ratio 
    //#include "ustar.H"
    //#include "1dProfileNut.H"
    //#include "1dProfileWind.H"
    //#include "computeUstarThetastar.H" // Compute friction coefficients 
    //#include "1dProfileTemp.H" // One-dimensional vertical Temperature profile
#include "1dSolver.H"
    //#include "1dProfileqV.H" // One-dimensional mixing ratio profile 
#include "createWindandTurbulence.H" // Create the three-dimensional fields 
#include "createTemperatureandHumidity.H"    
#include "adjustTurbulenceConstants.H"
    scalar radiationcounter=0;
    while (pimple.loop())
      {
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
	      //	      if(radiationcounter == callRadiation)
	      if(radiationcounter==15 || radiationcounter==30 || radiationcounter==45 
		 || radiationcounter==90 || radiationcounter==120 || radiationcounter==150 || radiationcounter==210)
		{
		  radiation->correct();
		  //		  radiationcounter=0;
		  firstTime=false;
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
	volScalarField diffuseSW=0*T;
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
	/*if(hotSpot && fastCounter>50)
	  {
#include "UEqn.H"
#include"nutEqn.H"
#include "pEqn.H"
} */
	fastCounter=fastCounter+1;
      }
    {
    // PMV
      //      #include "fvDOM.H"
      //    const radiationModel& MRTradiation =db().lookupObject<radiationModel>("radiationProperties");
    //    const fvDOM& dom(refCast<const fvDOM>(MRTradiation));
    }
    runTime.write();
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
	<< nl << endl;
  }
	

