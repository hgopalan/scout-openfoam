//scalar TZRef_(ABLDict.lookupOrDefault<scalar>("TZRef",2.0));
scalar ZREF_(ABLDict.lookupOrDefault<scalar>("ZRef",10.0));
scalar Z0_(ABLDict.lookupOrDefault<scalar>("z0",1e-4));
scalar mAvail_(ABLDict.lookupOrDefault<scalar>("mAvail",0.1));
Switch Nudge(ABLDict.lookupOrDefault<Switch>("Nudge",false));
Switch Damp(ABLDict.lookupOrDefault<Switch>("Damp",true));
scalar callRadiation(ABLDict.lookupOrDefault<scalar>("callRadiation",50));
//scalar timeScale_(ABLDict.lookupOrDefault<scalar>("timeScale",3600));
Switch writeSummary(ABLDict.lookupOrDefault<Switch>("writeSummary",true));
word solverType_(ABLDict.lookupOrDefault<word>("solverType","steadyHumid"));
//scalar counter=0;
scalar tempcounter=0;
scalar kappa_=0.41;
scalar locU1=interpolateXY(U.time().value(),timevalues,uvalues);
scalar locV1=interpolateXY(U.time().value(),timevalues,vvalues);
scalar M=Foam::sqrt(sqr(locU1)+sqr(locV1));
vector locFlowDir(locU1/M,locV1/M,0);
//scalar surfacez0=Z0_;
scalar MOL=interpolateXY(U.time().value(),timevalues,movalues);
if(Foam::mag(MOL)>25)
  MOL=MOL;
 else if(Foam::neg(MOL) && MOL> -25)
   MOL=-25;
 else if(Foam::pos(MOL) && Foam::mag(MOL)<25)
   MOL=25;

scalar PBLH=interpolateXY(U.time().value(),timevalues,pblhvalues);
scalar MBL=interpolateXY(U.time().value(),timevalues,mblvalues);
scalar RH=interpolateXY(mesh.time().value(),timevalues,rhvalues);
TRef=interpolateXY(U.time().value(),timevalues,TRefvalues);
scalar UUpper=interpolateXY(mesh.time().value(),timevalues,UUppervalues);
scalar VUpper=interpolateXY(mesh.time().value(),timevalues,VUppervalues);
scalar TUpper=min(interpolateXY(mesh.time().value(),timevalues,TUppervalues),306);
scalar qVUpper=interpolateXY(mesh.time().value(),timevalues,qVUppervalues);
solverType solverTypeModel_(solverTypeNames_.get("solverType", ABLDict));
ZREF_=ZREF_+Z0_;
Info<<"Initializing Flow Solver"<<endl;
