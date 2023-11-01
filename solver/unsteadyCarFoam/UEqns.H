surfaceScalarField alphaRain=fvc::interpolate(alpharain[phase_no]);
volTensorField gradU=fvc::grad(Urain[phase_no]);
volScalarField sourceU= treeDragCoeff*mag(Urain[phase_no]);
fvVectorMatrix UrainEqn
(
 fvm::div(alphaRain*phirain[phase_no],Urain[phase_no],"div(phiRain,Urain)")
==
 fvm::laplacian(alpharain[phase_no]*nutrain,Urain[phase_no])
 +fvc::div(alpharain[phase_no]*nutrain*dev2(gradU.T()),"div((alpha*(-nutrain*grad(U).T())))")
 -fvc::div(alpharain[phase_no]*(2.0/3.0)*I*nutrain*S/Foam::sqrt(0.09),"div((alpha*(-nutrain*grad(U).T())))")
 +alpharain[phase_no]*g
 + alpharain[phase_no]*((3* mua * CdRe)/(4 * rhop * dp * dp))*UCheat
 - fvm::Sp(((3* mua * CdRe)/(4 * rhop * dp * dp))*alpharain[phase_no],Urain[phase_no])
 -fvm::Sp(sourceU,Urain[phase_no])
 );
UrainEqn.relax();
UrainEqn.solve();
#include "URainConstraint.H"
Urain[phase_no].correctBoundaryConditions();

//scalar counter=0;

phirain[phase_no] = fvc::interpolate(Urain[phase_no]) & mesh.Sf();
//Urain[phase_no].correctBoundaryConditions();

/*Info<<"Failed cells count:"<<counter<<endl; */
Info<<"Rain Velocity: "<<gMax(Urain[phase_no])<<"    "<<gMin(Urain[phase_no])<<endl;
Info<<"Wind Velocity: "<<gMax(U)<<"   "<<gMin(U)<<endl;
