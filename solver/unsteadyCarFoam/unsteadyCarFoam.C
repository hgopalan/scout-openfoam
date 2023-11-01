/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
  \\/     M anipulation  |
  -------------------------------------------------------------------------------
  License
  This file is part of OpenFOAM.

  OpenFOAM is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
*/
// HG - Update comments sections 



#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "fvc.H"
#include "fvm.H"
//#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "scalarIOList.H"
#include "wallFvPatch.H"
#include "zeroGradientFvPatchField.H"
#include "interpolateXY.H"
#include "bound.H"
#include "wallDist.H"
#include "coordinateSystem.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "IFstream.H"
#include "OFstream.H"
#include "updateValues_V2.H"



int main(int argc, char *argv[])
{


    #include "defParameters.H"
  argList args(argc, argv);  

#include "postProcess.H"
  //#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createControl.H"
#include "createFields.H"
#include "createFvOptions.H"
#include "createTimeControls.H"
#include "CourantNo.H"
#include "setInitialDeltaT.H"
#include "initContinuityErrs.H"





  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // Initialize
#include "initializeSolver.H" // Lot of variables
#include "1dInitalize.H" // Create vertical scalar list for Temperature and mixing ratio 
  //#include "createWindandTurbulence.H" // Create the three-dimensional fields 
  //#include "computeUstarThetastar.H" // Compute friction coefficients 
  //#include "1dProfileTemp.H" // One-dimensional vertical Temperature profile
  //#include "1dProfileqV.H" // One-dimensional mixing ratio profile 
    //#include "createTemperatureandHumidity.H"
#include "pmvinitialize.H" // Initialize the three-dimensional fields for thermal comfort 
  // Setup Nudge Domain 
  if(Nudge)
#include "markNudgeCells.H"
volVectorField nudgeU_=0.0*U;
volScalarField nudgeNut_=0.0*nut;    
    scalar nu=1e-5;
    scalar turnOffShading=1;
    // Option 1 - Steady No Temperature
    switch (solverTypeModel_)
      {
      case steadyWind:
	{
	  turnOffShading=0;
#include "runWindSolver.H"
	  break;
	}
      case unsteadyWind:
	{
	  turnOffShading=0;
#include "runWindSolver.H"
	  break;
	  } 
      case steadyThermal:
	{
#include "runThermalSolver.H"
	  break;
	}
      case unsteadyThermal:
	{
#include "runThermalSolver.H"
	  break;
	  } 
      case steadyHumid:
	{
#include "runHumidSolver.H"
	  break;
	}
      case unsteadyHumid:
	{
#include "runHumidSolver.H"
	  break;
	}
      }
  Info<< "End\n" << endl;
  return 0;
}



// ************************************************************************* //
