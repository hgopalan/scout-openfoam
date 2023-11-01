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
