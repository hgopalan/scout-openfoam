Info<< nl << "Calculating wet spots\n" << endl;

volScalarField wet
(
 IOobject
 (
  "totalWetPercent",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE		
  ),
 mesh,
 dimensionedScalar("wet",dimensionSet(0,0,0,0,0,0,0),0)
 );
for (label i=0;i<phases.size();i++)
  {
    string phase_scrname = "diameter_";
    string phase_no_str;         
    std::stringstream convert;
    convert << i+1;
    //convert<<defaultPhases[i]*1000;
    phase_no_str = convert.str();
    phase_scrname.append(phase_no_str);
    //    phase_scrname.append("_m");
    /*volScalarField alphatemp
      (
      IOobject
      (
      "wet",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
      alpharain[i]
      );
	
      scalar maxalpha = max(alphatemp).value();
	
      forAll(alphatemp,celli) 
      {
      if(alphatemp[celli] > scalar(0.0001*maxalpha))
      {
      wet[celli] = wet[celli]+scalar(1);
      wet[celli]=min(wet[celli],1);
      }
      }	
      for (label j=0;j<wet.boundaryField().size();j++)
      {
      forAll(alphatemp.boundaryField()[j],celli)
      if (alphatemp.boundaryField()[j][celli] > scalar(0.0001*maxalpha))
      {
      wet.boundaryFieldRef()[j][celli] = wet.boundaryFieldRef()[j][celli]+scalar(1);
      wet.boundaryFieldRef()[j][celli] = min(wet.boundaryFieldRef()[j][celli],1);
      }

      } */
    scr.set
      (
       i,	
       new volScalarField 
       (
	IOobject
	(
	 phase_scrname,
	 runTime.timeName(),
	 mesh,
	 IOobject::NO_READ,
	 IOobject::NO_WRITE
	 ),
	mesh,
	dimensionedScalar("wet",dimensionSet(0,0,0,0,0,0,0),1)
	)
       );
    scr[i]=alpharain[i]/gMax(alpharain[i])*100;
    if(i==0)
      wet=alpharain[i]/gMax(alpharain[i])*100;
    else
      wet=max(wet,alpharain[i]/gMax(alpharain[i])*100);
    /*forAll(alphatemp,celli) 
      {
      if(alphatemp[celli] < scalar(0.0001*maxalpha))
      scr[i][celli] = scalar(0);
      }
      for (label j=0;j<wet.boundaryField().size();j++)
      {
      forAll(alphatemp.boundaryField()[j],celli)
      if (alphatemp.boundaryField()[j][celli] < scalar(0.0001*maxalpha))
      scr[i].boundaryFieldRef()[j][celli] = scalar(0);
      }*/
    scr[i].write();
  }
//wet=wet;///phases.size();
Info << nl << "Done calculating wet cells " << endl;
Info << nl << max(wet) << " " << min(wet) << endl;
wet.write();

