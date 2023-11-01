
/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
  \\/     M anipulation  |
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

#ifndef ABLMOWallTFvPatchScalarField_H
#define ABLMOWallTFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
  !----------------------------------------------------------------------------!
  ! Surface energy Balance for Temperature                        !
  !----------------------------------------------------------------------------!                                                    
  This routine implements the surface energy balance for temperature. The 
  surface energy balance is given by (neglecting heat storage)
  Short-Wave Radiation + Long-Wave Radiation = Sensible + Latent + Ground.
  Currently latent heat-flux is neglected. Long-wave radiation has temperature 
  dependence and the equation is non-linear. To avoid this, the temperature 
  used for the long-wave radiation is lagged. Also the boundary condition is 
  implemented as a gradient boundary condition for roubstness. As a result
  ground heat-flux is also lagged. The final form is 

  SW + LW - G = rho Cp nut/Prt dT/dn
  
  DT/dn is fixed on the walls. When the friction velocity ustar ->0, nut can 
  become zero. To avoid this issue, the value of nut is modified. See the .C 
  file for more information. 
 

  References: 

  Log File:
  1. 28/3/2019: Adopted for OF-v1812. 

  Authors:
  Harish Gopalan (gopalanh@ihpc.a-star.edu.sg)
*/ 

/* Usage 
   terrain
   {
   type ABLMOWallT;
   z0 uniform 1.000000; // Roughness length in m. Default is 0.1 m. 
   value uniform 300.15; // This is a placeholder. Use values in K. 
   }
*/

namespace Foam
{

  /*---------------------------------------------------------------------------*\
    Class ABLMOWallTFvPatchScalarField Declaration
    \*---------------------------------------------------------------------------*/

  class ABLMOWallTFvPatchScalarField
    :
    public fixedValueFvPatchScalarField
  {
  protected:

    // Protected data
    //- Surface roughness length
    // Heat-flux methods: (1) viewFactor (2) solarLoad (3) constant
    // Currently turned off due to some confusing decomposePar issues
    // We assume that the long-wave method is fvDOM. If not long-wave 
    // contribution is neglected. 
    //      word surfaceHeatFlux_; // Reading heat-flux method 
    scalar z0_; // Roughness length in m
    scalarList matDensity_; // Density of each layer in kg/m^3
    scalarList specificHeat_; // Heat capacity of each layer 
    scalarList thermalCon_; // Thermal Conductivity of each layer in W/m^2
    scalarList thickness_; // Thickness of each layer in m
    scalarField T1_; // First layer
    scalarField T2_; // Second layer
    scalarField T3_; // Third layer
    scalarField T4_; // Fourth layer
    scalarField Tdeep_; // Soil or internal temperature 

    // Protected Member Functions


  public:

    //- Runtime type information
    TypeName("ABLMOWallT");


    // Constructors

    //- Construct from patch and internal field
    ABLMOWallTFvPatchScalarField
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&
     );

    //- Construct from patch, internal field and dictionary
    ABLMOWallTFvPatchScalarField
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const dictionary&
     );

    //- Construct by mapping given
    //  ABLMOWallTFvPatchScalarField
    //  onto a new patch
    ABLMOWallTFvPatchScalarField
    (
     const ABLMOWallTFvPatchScalarField&,
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const fvPatchFieldMapper&
     );

    //- Construct as copy
    ABLMOWallTFvPatchScalarField
    (
     const ABLMOWallTFvPatchScalarField&
     );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
      return tmp<fvPatchScalarField>
	(
	 new ABLMOWallTFvPatchScalarField(*this)
	 );
    }

    //- Construct as copy setting internal field reference
    ABLMOWallTFvPatchScalarField
    (
     const ABLMOWallTFvPatchScalarField&,
     const DimensionedField<scalar, volMesh>&
     );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone
    (
     const DimensionedField<scalar, volMesh>& iF
     ) const
    {
      return tmp<fvPatchScalarField>
	(
	 new ABLMOWallTFvPatchScalarField(*this, iF)
	 );
    }


    // Member functions

    // Acces functions

    // Return z0


    // Mapping functions

    //- Map (and resize as needed) from self given a mapping object
    virtual void autoMap(const fvPatchFieldMapper&);

    //- Reverse map the given fvPatchField onto this fvPatchField
    virtual void rmap
    (
     const fvPatchScalarField&,
     const labelList&
     );

    tmp<scalarField> phiM(scalarField &p,scalar& MOL_) const;
    tmp<scalarField> phiH(scalarField &p,scalar& MOL_) const;
    tmp<scalarField> psiM(scalarField &p,scalar& MOL_) const;
    tmp<scalarField> psiH(scalarField &p,scalar& MOL_) const;
    tmp<scalarField> totalRadiation() const;
    tmp<scalarField> latentHeat(scalarField& wstar,scalarField& nhat,scalar MOL) const;
    tmp<scalarField> groundHeat() const;
    tmp<scalarField> firstLayer(scalarField &Twall,scalar dt) const;
    tmp<scalarField> secondLayer(scalar dt) const;
    tmp<scalarField> thirdLayer(scalar dt) const;
    tmp<scalarField> fourthLayer(scalar dt) const;
    void minMax(scalarField& field) const;
    // I-O
    //- Update the coefficients associated with the patch field
    using fixedValueFvPatchScalarField::updateCoeffs;
    virtual void updateCoeffs();
    //- Write
    virtual void write(Ostream&) const;
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
