{
  // Evapotranspiration effects
  forAll(leafTemperature,cellI)
    {
      scalar rhoRef=rhok[cellI];
      scalar pRef=p_rgh[cellI]+rhok[cellI]*gh[cellI]+101325;
      //      scalar factor=0.622;
      scalar ra_=130*Foam::sqrt(0.5/(mag(U[cellI])+0.05));
      scalar rs_=300; //Temporary fix
      scalar pvsat=610.94*Foam::exp(17.625*(T[cellI]-273.15)/(T[cellI]-30.11));
      leafTemperature[cellI]=pos(treeblanking[cellI])*rhoRef/pRef*0.622/(ra_+rs_)*pvsat*(1-relHum[cellI]);
    }
  leafTemperature.correctBoundaryConditions();
  Info<<"Source:"<<gMin(leafTemperature)<<"  "<<gMax(leafTemperature)<<endl;
}
