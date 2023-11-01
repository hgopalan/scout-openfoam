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
    phase_no_str = convert.str();
    phase_scrname.append(phase_no_str);
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
    //    if(i==0)
    //  wet=alpharain[i]/gMax(alpharain[i])*100;
    //else
    //  wet=max(wet,alpharain[i]/gMax(alpharain[i])*100);
    scr[i].write();
  }
//wet=wet;///phases.size();
scalar gmax=0;
for (label i=0;i<phases.size();i++)
  {
    //gmax=gmax+gMax(alpharain[i]);
    wet=wet+alpharain[i]/gMax(alpharain[i])*100;
  }
wet=wet/phases.size();
Info << nl << "Done calculating wet cells " << endl;
Info << nl << max(wet) << " " << min(wet) << endl;
wet.write();

