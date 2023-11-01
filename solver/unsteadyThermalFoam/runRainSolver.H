scalar raincounter=0;
while(raincounter<=250)
  {
    Info<<"Rain Iteration number:"<<raincounter<<endl;
    for (label phase_no = 0; phase_no < phases.size(); phase_no++)
      {
	if(raincounter==0)
	  {
	    forAll(Urain[phase_no],cellI)
	  {
	    Urain[phase_no][cellI][0]=UCheat[cellI][0];
	    Urain[phase_no][cellI][1]=UCheat[cellI][1];
	    Urain[phase_no][cellI][2]=-Vt[phase_no];
	  }
	  Urain[phase_no].correctBoundaryConditions();
	  phirain[phase_no]=fvc::interpolate(Urain[phase_no]) & mesh.Sf();
      }
       
    // phi is used by the inletOutlet boundary condition and courant number calculation
    surfaceScalarField phi
      (
       IOobject
       (
	"phi",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
	),
       phirain[phase_no]
       );
#include "alphaEqns.H"			
    dimensionedScalar dp ("dp", dimensionSet(0,1,0,0,0,0,0), phases[phase_no]);
    volScalarField magUr = mag(UCheat - Urain[phase_no]);			
    Re = (magUr*dp*rhoa)/mua;
    CdRe = GET_CdRe(Re);
    CdRe.correctBoundaryConditions();
    //Info<<"Solving turbulence"<<endl;
    volScalarField S=(sqrt(2*(symm(fvc::grad(U)) && fvc::grad(U)))+1e-15);
    tfl = 0.2/(Foam::sqrt(0.09)*S+SMALL);
    tp = (4*rhop*dp*dp)/(3*mua*CdRe+SMALL);
    Ctrain[phase_no] = (tfl/(tfl+tp));				
    nutrain = nut*Foam::sqr(Ctrain[phase_no]);
#include "UEqns.H"
    Info<<" Rain Turbulence:"<<gMax(nutrain)<<"  "<<gMin(nutrain)<<endl;
    Info<<" Wind Turbulence:"<<gMax(nut)<<"  "<<gMin(nut)<<endl;
  }
raincounter=raincounter+1;
}
