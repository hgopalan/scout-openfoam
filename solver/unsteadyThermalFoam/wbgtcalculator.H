/*const volScalarField swr_=U.db().lookupObject<volScalarField>("qmrtsw");
scalar S=directSolarRad_+diffuseSolarRad_;
scalar fdir=directSolarRad_/(S+1E-5);
scalar fdif=diffuseSolarRad_/(S+1E-5);
volScalarField shading=swr_/(gMax(swr_)+1e-5);
shading=max(shading,0*shading);
volScalarField Td=T-(100-RH)/5-273.15;//Crude
volScalarField Ta=T-273.15;
volScalarField ea=Foam::exp(17.67*(Td-Ta)/(Td+243.5))*(1.0007+0.00000346*(p+101325))*Foam::exp(17.5*Ta/(Ta+240.97));
volScalarField emissA=0.58*Foam::pow(ea,1.0/7.0);
volScalarField wbgt_b=S*shading*(fdir/(2.268e-8*Foam::cos(beta_)+1e-5)+1.2*fdif/5.67e-8)+emissA*Foam::pow(Ta,4);
volScalarField wind=mag(U)*3600;
wind=max(wind,1+0*wind);
volScalarField wbgt_c=0.315*Foam::pow(wind,0.58)/5.3865e-8;
volScalarField Tg=(wbgt_b+wbgt_c*(T-273.15)+7680000)/(wbgt_c+256000);
volScalarField Twet=Ta*Foam::atan(0.16*Foam::sqrt(relHum*100+8.313))+Foam::atan(Ta+relHum*100)-Foam::atan(relHum*100-1.676)+0.0039*Foam::pow(relHum*100,1.5)*Foam::atan(0.023*(relHum*100))-4.686;
WBGT=0.7*Twet+0.2*Tg+0.1*Ta;
WBGT.correctBoundaryConditions();
Info<<"WBGT:"<<gMin(WBGT)<<" "<<gMax(WBGT)<<endl;
*/
// Tested for Singapore https://www.wbgt.sg/#chart
forAll(WBGT,cellI)
{
  if(relHum[cellI]<=0.4)
    WBGT[cellI]=T[cellI]-273.15;
  else
    WBGT[cellI]=(relHum[cellI]-0.40)/0.6*10+(T[cellI]-273.15-2);
}
WBGT.correctBoundaryConditions();
Info<<"WBGT:"<<gMin(WBGT)<<" "<<gMax(WBGT)<<endl;
