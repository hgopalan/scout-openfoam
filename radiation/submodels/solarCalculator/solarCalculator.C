/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solarCalculator.H"
#include "Time.H"
#include "unitConversion.H"
#include "constants.H"
#include "interpolateXY.H"
using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solarCalculator, 0);
}


const Foam::Enum
<
    Foam::solarCalculator::sunDirModel
>
Foam::solarCalculator::sunDirectionModelTypeNames_
({
    { sunDirModel::mSunDirConstant, "sunDirConstant" },
    { sunDirModel::mSunDirTracking, "sunDirTracking" },
    { sunDirModel::mSunDirTrackingAdv, "sunDirTrackingAdv" },
});


const Foam::Enum
<
    Foam::solarCalculator::sunLModel
>
Foam::solarCalculator::sunLoadModelTypeNames_
({
    { sunLModel::mSunLoadConstant, "sunLoadConstant" },
    {
        sunLModel::mSunLoadFairWeatherConditions,
        "sunLoadFairWeatherConditions"
    },
    { sunLModel::mSunLoadTheoreticalMaximum, "sunLoadTheoreticalMaximum" },
    { sunLModel::mLiuJordan, "LiuJordan" },     
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solarCalculator::calculateBetaTheta()
{
    scalar runTime = 0.0;
    switch (sunDirectionModel_)
    {
        case mSunDirTracking:
        {
            runTime = mesh_.time().value();
            break;
        }
        case mSunDirConstant:
        {
            break;
        }
        case mSunDirTrackingAdv:
        {
            runTime = mesh_.time().value();
            break;
        }
    }

    scalar LSM = 15.0*(dict_.get<scalar>("localStandardMeridian"));

    scalar D = dict_.get<scalar>("startDay") + runTime/86400.0;
    scalar M = 6.24004 + 0.0172*D;
    scalar EOT = -7.659*sin(M) + 9.863*sin(2*M + 3.5932);

    dict_.readEntry("startTime", startTime_);

    scalar LST =  startTime_ + runTime/3600.0;

    scalar LON = dict_.get<scalar>("longitude");

    scalar AST = LST + EOT/60.0 + (LON - LSM)/15;

    scalar delta = 23.45*sin(degToRad((360*(284 + D))/365));

    scalar H = degToRad(15*(AST - 12));

    scalar L = degToRad(dict_.get<scalar>("latitude"));

    scalar deltaRad = degToRad(delta);
    beta_ = max(asin(cos(L)*cos(deltaRad)*cos(H) + sin(L)*sin(deltaRad)), 1e-3);
    theta_ = acos((sin(beta_)*sin(L) - sin(deltaRad))/(cos(beta_)*cos(L)));

    // theta is the angle between the SOUTH axis and the Sun
    // If the hour angle is lower than zero (morning) the Sun is positioned
    // on the East side.
    if (H < 0)
    {
        theta_ += 2*(constant::mathematical::pi - theta_);
    }

    if (debug)
    {
        Info << tab << "altitude : " << radToDeg(beta_) << endl;
        Info << tab << "azimuth  : " << radToDeg(theta_) << endl;
    }
}


void Foam::solarCalculator::calculateSunDirection()
{
    gridUp_  = normalised(dict_.get<vector>("gridUp"));
    eastDir_ = normalised(dict_.get<vector>("gridEast"));

    coord_.reset
    (
        new coordinateSystem("grid", Zero, gridUp_, eastDir_)
    );

    // Assuming 'z' vertical, 'y' North and 'x' East
    direction_.z() = -sin(beta_);
    direction_.y() =  cos(beta_)*cos(theta_); // South axis
    direction_.x() =  cos(beta_)*sin(theta_); // West axis

    direction_.normalise();

    if (debug)
    {
        Info<< "Sun direction in absolute coordinates : " << direction_ <<endl;
    }

    // Transform to actual coordinate system
    direction_ = coord_->transform(direction_);

    if (debug)
    {
        Info<< "Sun direction in the Grid coordinates : " << direction_ <<endl;
    }
}


void Foam::solarCalculator::init()
{
    switch (sunDirectionModel_)
    {
        case mSunDirConstant:
        {
            if (dict_.readIfPresent("sunDirection", direction_))
            {
                direction_.normalise();
            }
            else
            {
                calculateBetaTheta();
                calculateSunDirection();
            }

            break;
        }
        case mSunDirTracking:
        {
            if (word(mesh_.ddtScheme("default")) == "steadyState")
            {
                FatalErrorInFunction
                    << " Sun direction model can not be sunDirtracking if the "
                    << " case is steady " << nl << exit(FatalError);
            }

            dict_.readEntry
            (
                "sunTrackingUpdateInterval",
                sunTrackingUpdateInterval_
            );

            calculateBetaTheta();
            calculateSunDirection();
            break;
        }
        case mSunDirTrackingAdv:
        {
            if (word(mesh_.ddtScheme("default")) != "steadyState")
            {
             dict_.readEntry
             (
                 "sunTrackingUpdateInterval",
                 sunTrackingUpdateInterval_
             );
            }
            calculateBetaTheta();
            calculateSunDirection();
            directSolarRad_=interpolateXY(mesh_.time().value(),timeprofile_,directswprofile_);
            diffuseSolarRad_=interpolateXY(mesh_.time().value(),timeprofile_,diffuseswprofile_);
            Info<<"Initial Radiation:"<<mesh_.time().value()<<"  "<<directSolarRad_<<"   "<<diffuseSolarRad_<<endl;            
            break;
		}
    }
    switch (sunLoadModel_)
    {
        case mSunLoadConstant:
        {
            dict_.readEntry("directSolarRad", directSolarRad_);
            dict_.readEntry("diffuseSolarRad", diffuseSolarRad_);
            break;
        }
        case mSunLoadFairWeatherConditions:
        {
            dict_.readIfPresent
            (
                "skyCloudCoverFraction",
                skyCloudCoverFraction_
            );

            dict_.readEntry("A", A_);
            dict_.readEntry("B", B_);

            if (!dict_.readIfPresent("beta", beta_))
            {
                calculateBetaTheta();
            }

            directSolarRad_ =
                (1.0 - 0.75*pow(skyCloudCoverFraction_, 3.0))
              * A_/exp(B_/sin(beta_));

            dict_.readEntry("groundReflectivity", groundReflectivity_);
            break;
        }
        case mSunLoadTheoreticalMaximum:
        {
            dict_.readEntry("Setrn", Setrn_);
            dict_.readEntry("SunPrime", SunPrime_);
            directSolarRad_ = Setrn_*SunPrime_;

            dict_.readEntry("groundReflectivity", groundReflectivity_);
            break;
        }
        case mLiuJordan:
        {
			break;
		}
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarCalculator::solarCalculator
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dict),
    direction_(Zero),
    directSolarRad_(0.0),
    diffuseSolarRad_(0.0),
    groundReflectivity_(0.0),
    A_(0.0),
    B_(0.0),
    beta_(0.0),
    theta_(0.0),
    skyCloudCoverFraction_(0.0),
    Setrn_(0.0),
    SunPrime_(0.0),
    C_(dict.get<scalar>("C")),
    sunDirectionModel_
    (
        sunDirectionModelTypeNames_.get("sunDirectionModel", dict)
    ),
    sunLoadModel_(sunLoadModelTypeNames_.get("sunLoadModel", dict)),
    coord_(),
  profilesize(dict.lookupOrDefault<scalar>("profilesize",1000)),
  firsttime_(true),
  timeprofile_(profilesize,-1e30),
  directswprofile_(profilesize,0.0),
  diffuseswprofile_(profilesize,0.0) 
{
  if(firsttime_)
    {
      firsttime_=false;
      List <List<scalar> > profileTable(dict.lookup("profileTable"));
      scalarField zProfile(profileTable.size(),0.0);
      forAll(zProfile,i)
	{
	  timeprofile_[i] = profileTable[i][0];
	  directswprofile_[i] = profileTable[i][1];
	  diffuseswprofile_[i] = profileTable[i][2];
	}
    }
    init();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solarCalculator::correctSunDirection()
{
    switch (sunDirectionModel_)
    {
        case mSunDirConstant:
        {
            break;
        }
        case mSunDirTracking:
        {
            calculateBetaTheta();
            calculateSunDirection();
            directSolarRad_ = A_/exp(B_/sin(max(beta_, ROOTVSMALL)));
            break;
        }
        case mSunDirTrackingAdv:
        {
           if (word(mesh_.ddtScheme("default")) == "steadyState")
            {
			 break;
            }
            else
            {
            calculateBetaTheta();
            calculateSunDirection();
            directSolarRad_=interpolateXY(mesh_.time().value(),timeprofile_,directswprofile_);
            diffuseSolarRad_=interpolateXY(mesh_.time().value(),timeprofile_,diffuseswprofile_);
            Info<<"Initial Radiation:"<<mesh_.time().value()<<"  "<<directSolarRad_<<"   "<<diffuseSolarRad_<<endl;            
            break;
		}
		}        
    }
}


// ************************************************************************* //
