scalar LSM = 15.0*(solarDict.get<scalar>("localStandardMeridian"));
scalar D = solarDict.get<scalar>("startDay") + U.time().value()/86400.0;
scalar M = 6.24004 + 0.0172*D;
scalar EOT = -7.659*Foam::sin(M) + 9.863*Foam::sin(2*M + 3.5932);
scalar startTime_=solarDict.get<scalar>("startTime");
scalar LST =  startTime_ + U.time().value()/3600.0;
scalar LON = solarDict.get<scalar>("longitude");
scalar AST = LST + EOT/60.0 + (LON - LSM)/15;
scalar delta = 23.45*Foam::sin(((360*(284 + D))/365)*3.14/180);
scalar H = (15*(AST - 12))*3.14/180;
scalar L = (solarDict.get<scalar>("latitude"))*3.14/180;
scalar deltaRad = delta*3.14/180;
scalar beta_ = max(Foam::asin(Foam::cos(L)*Foam::cos(deltaRad)*Foam::cos(H) + Foam::sin(L)*Foam::sin(deltaRad)), 1e-3);
scalar theta_ = Foam::acos((Foam::sin(beta_)*Foam::sin(L) - Foam::sin(deltaRad))/(Foam::cos(beta_)*Foam::cos(L)));
vector gridUp_  = normalised(solarDict.get<vector>("gridUp"));
vector eastDir_ = normalised(solarDict.get<vector>("gridEast"));
vector sunDirection_;
autoPtr<coordinateSystem> coord_;
coord_.reset (new coordinateSystem("grid", Zero, gridUp_, eastDir_));
// Assuming 'z' vertical, 'y' North and 'x' East
sunDirection_.z() = - Foam::sin(beta_);
sunDirection_.y() =  Foam::cos(beta_)*Foam::cos(theta_); // South axis
sunDirection_.x() =  Foam::cos(beta_)*Foam::sin(theta_); // West axis
sunDirection_.normalise();
sunDirection_ = coord_->transform(sunDirection_);
scalar directSolarRad_=interpolateXY(U.time().value(),radTimeProfile_,radDirectswProfile_);
scalar diffuseSolarRad_=interpolateXY(U.time().value(),radTimeProfile_,radDiffuseswProfile_);
vector verticalDir(0,0,1);
forAll(diffuseSW,cellI)
{
  if(pos(treeblanking[cellI]))
    {
      vector cc=mesh.C()[cellI];
      scalar eta=cc[2]/max(cc[2],treeShading[cellI]);
      diffuseSW[cellI]=diffuseSolarRad_*Foam::exp(-10*(1-eta));
    }
}
Info<<" Sun is here:"<<sunDirection_<<endl;
#include "wbgtcalculator.H"
