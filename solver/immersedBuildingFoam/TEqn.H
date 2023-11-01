{
#include "readSolarCalculator.H" // Computes Diffuse solar load penetration on trees 
  tempcounter=0;
  //  Info<<"Before AlphaT"<<endl;
  alphaT = nut/Prt+buildingBlanking*nu;
  alphaT.correctBoundaryConditions();
  //  Info<<"After AlphaT"<<endl;
  volScalarField alphaEff("alphaEff",alphaT);
    // Long-Wave (sky)
  //  volScalarField cosEpsilon=(sunDirection_ & -mesh.C()/mag(mesh.C()));
  // Incoming Radiation from Clouds (not accounted in the exchanges)
  //  volScalarField qSky_=0.5*(1+cosEpsilon)*(213+5.5*(T-273.15));
  scalar Tsky_=0.0555*Foam::pow(TRef,1.5);
  //  volScalarField qSky_=0.5*(1+cosEpsilon)*5.67e-8*(Foam::pow(Tsky_,4)-Foam::pow(T,4));
  volScalarField qSky_=5.67e-8*(Foam::pow(Tsky_,4)-Foam::pow(T,4));
  // Diffuse Radiation from Sun 
  // Convert to Potential Temperature
  surfaceScalarField phiGhost=0*phi;
  forAll(buildingBlankingSurface,faceI)
  {
	  if(buildingBlankingSurface[faceI]<0.2 && abs(phi[faceI])>0)
		phiGhost[faceI]=-phi[faceI];
	}
  fvScalarMatrix TEqn
    (
     fvm::ddt(T)
     ==
      -fvm::div(phi-phiGhost,T)
      +fvm::laplacian(alphaEff,T)
     +neg(treeblanking)*radiation->ST(rhok*cpRef, T)
     +treeLAD*pos(treeblanking)*(radiation->ST(rhok*cpRef, T)+(diffuseSW+qSky_)/(rhok*cpRef)-2.2e6/1000.0*leafTemperature)
     +pos(sourceBlanking)*(temperatureSource/(rhok*cpRef))
     -(1-buildingBlanking)*fvm::laplacian(2.25/1e6+0*alphaEff,T)
     ); 
  TEqn.relax();
     //+fvOptions(T)  
  //fvOptions.constrain(TEqn);
  TEqn.setValues(blankCellList,internalTemperature);  
  TEqn.solve();
  // In the unphysical gaps in buildings temperature can built-up.
  // Limiting avoids the built up. 
  scalar lowerT=TRef-8;
  scalar upperT=TRef+8;
  forAll(T,cellI)
    {
      if(T[cellI]<lowerT)
	{
	  tempcounter=tempcounter+1;
	  T[cellI]=lowerT;
          scalar esat=610.094*Foam::exp(17.625*(T[cellI]-273.15)/(T[cellI]-30.11));
          scalar qsat=(0.622 *esat/(101325-0.378*esat));
          qv[cellI]=relHum[cellI]*qsat;
	}
      if(T[cellI]>upperT)
	{
          tempcounter=tempcounter+1;
          T[cellI]=pos(unphysicalCells[cellI]-3)*upperT+neg(unphysicalCells[cellI]-3)*TRef;
	  scalar esat=610.094*Foam::exp(17.625*(T[cellI]-273.15)/(T[cellI]-30.11));
	  scalar qsat=(0.622 *esat/(101325-0.378*esat));
	  qv[cellI]=relHum[cellI]*qsat;
	}
    }
  reduce(tempcounter, sumOp<scalar>());
  if(tempcounter>0)
    {
      T.correctBoundaryConditions();
      qv.correctBoundaryConditions();
    }
  //  if(pimple.finalIter())                                                                     
  //    radiation->correct(); 
  //  radiation->correct();
  // TV
  thetaV=T*(1+(461.51/287-1)*qv);
  //  scalar cellesat=610.94*Foam::exp(17.625*(TRef-273.15)/(TRef-30.11));
  //  scalar cellqsat=(0.622 *cellesat/(101325-0.378*cellesat));
  //  scalar qvRef=RH*cellqsat;
  //  TRef=TRef*(1+(461.51/287-1)*qvRef);
  rhok= 1 -buildingBlanking/volTRef*(T-volTRef);
  TRef=interpolateXY(mesh.time().value(),timevalues,TRefvalues);
  Info<<"Potential Temperature Bounds:"<<gMin(T)<<"    "<<gMax(T)<<endl;                                             
  Info<<"Virtual Potential Temperature Bounds:"<<gMin(thetaV)<<"    "<<gMax(thetaV)<<endl;                           
  Info<<"Out of bound cell count for T:"<<tempcounter<<"  "<<TRef<<endl;      
}
 
