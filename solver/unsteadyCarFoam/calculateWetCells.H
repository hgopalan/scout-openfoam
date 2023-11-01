Info<< nl << "Calculating wet spots\n" << endl;

volScalarField wet
(
	IOobject
	(
		"wet",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("wet",dimensionSet(0,0,0,0,0,0,0),1)
);


for (int phase_no = 1; phase_no < phases.size(); phase_no++)
{
	volScalarField alphatemp
	(
		IOobject
		(
			"wet",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		alpharain[phase_no]
	);
	
	double maxalpha = max(alphatemp).value();
	
	forAll (alphatemp,celli) {
		if(alphatemp[celli] < scalar(0.0001*maxalpha)){
			wet[celli] = scalar(0);
		}
	}

	for (int i=0;i<wet.boundaryField().size();i++){
		forAll (alphatemp.boundaryField()[i],celli){
			if (alphatemp.boundaryField()[i][celli] < scalar(0.0001*maxalpha))
				wet.boundaryField()[i][celli] = scalar(0);
		}
	}
}


Info << nl << "Done calculating wet cells " << endl;
Info << nl << max(wet) << " " << min(wet);
wet.write();
