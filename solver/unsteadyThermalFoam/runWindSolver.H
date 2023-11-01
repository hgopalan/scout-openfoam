Info<<" Starting Humid Solver"<<endl;
#include "adjustTurbulenceConstants.H"
turb_C3.write();
blankcell.write();
vector vertical(0,0,1.0);
Switch eulerianWDR=ABLDict.lookupOrDefault<bool>("eulerianWDR", false);
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
#include "1dSolver.H"
#include "createWindandTurbulence.H" // Create the three-dimensional fields 
      //#include "computeUstarThetastar.H" // Compute friction coefficients 
      //#include "1dProfileTemp.H" // One-dimensional vertical Temperature profile
      //#include "1dProfileqV.H" // One-dimensional mixing ratio profile 
#include "createTemperatureandHumidity.H"    
#include "adjustTurbulenceConstants.H"
    // Dry the air 
    vertqV=1e-5+0*vertqV;
    vertqVProfile_=1e-5+0*vertqVProfile_;
    qv=1e-5+0*qv;
    qv.correctBoundaryConditions();
    forAll(qv,cellI)
      {
        scalar esat=610.94*Foam::exp(17.625*(T[cellI]-273.15)/(T[cellI]-30.11));
        scalar qsat=(0.622*esat/(101325-0.378*esat));
        relHum[cellI]=qv[cellI]/qsat;
      }
    relHum.correctBoundaryConditions();
    Info<<"Specific Humidity Bounds:"<<gMin(qv)<<"  "<<gMax(qv)<<endl;
    Info<<"Relative  Humidity Bounds:"<<gMin(relHum)<<"  "<<gMax(relHum)<<endl;
    //    scalar radiationcounter=0;
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
	volScalarField diffuseSW=0*T;
#include "TEqn.H"
	//	#include "qvEqn.H"
	//#include "computeEvapotranspiration.H"
	//#include "pmvcalculator.H"
#include "UEqn.H"
#include"nutEqn.H"
#include "pEqn.H"
	// Call Eulerian Wind Driven Rain                                                                                  
	if(eulerianWDR && pimple.finalIter())
	  {
#include "callEulerianRain.H"
#include "runRainSolver.H"
	    if (runTime.outputTime())
	      {
		for (label phase_no = 0; phase_no < phases.size(); phase_no++)
		  {
		    Urain[phase_no].write();
		    alpharain[phase_no].write();
		  }
#include "calculateWet.H"
#include "calculateCatchRatio.H"
	      }

	  }
      }
    runTime.write();
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
	<< nl << endl;
  }
	

