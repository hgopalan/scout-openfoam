/*---------------------------------------------------------------------------*\
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
/*
  !----------------------------------------------------------------------------!
  ! Wind Inflow Conditions for Atmospheric Boundary Layer                      !
  !----------------------------------------------------------------------------!                                                         
  This module implements the log-law profile for specifying inflow wind speed 
  under all stratifications. The original formulation is based on (JEP) and the 
  OpenFOAM formulation is based on (HG2018). This BC can work with both steady
  and unsteady reference wind speeds. The reference values are not specified 
  in the BC but in a separate file "constant/ABLDict". 

  References: 
  (JEP) Pieterse, J. E.,"CFD Investigation of the Atmospheric Boundary Layer under 
  Different Thermal Stability Conditions",Masters Thesis,Stellenbosch University,2013
  (HG2018) Gopalan, H., "Evaluation of Wray-Agarwal turbulence model for simulation of 
  neutral and non-neutral atmospheric boundary layers",JWEIA, 2018: 182, 322-329.

  Log File:
  1. 28/3/2019: Adopted for OF-v1812. 

  Authors:
  Harish Gopalan (gopalanh@ihpc.a-star.edu.sg)
*/ 

/* Usage 
   inflow
   {
   type unsteadyABLInflowU;
   ZRef 15.000000; // Default value is 10.0 m 
   z0 1.000000; // Default value is 0.1 m 
   value uniform (0 0 0); // This is a placeholder. Reference wind speed 
   specified in constant/ABLDict file. 
   }
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef unsteadyABLInflowUFvPatchVectorField_H
#define unsteadyABLInflowUFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------*\
    Class unsteadyABLInflowUFvPatchVectorField Declaration
    \*---------------------------------------------------------------------------*/

  class unsteadyABLInflowUFvPatchVectorField
    :
    public fixedValueFvPatchVectorField

  {
    scalar kappa_; // von-Karman constant, 0.41 
    scalar ZRef_; // Reference height in m 
    scalar z0_; // Roughness length in m 
    scalar zd_; // Displacement height in m 
  public:

    //- Runtime type information
    TypeName("unsteadyABLInflowU");


    // Constructors

    //- Construct from patch and internal field
    unsteadyABLInflowUFvPatchVectorField
    (
     const fvPatch&,
     const DimensionedField<vector, volMesh>&
     );

    //- Construct from patch, internal field and dictionary
    unsteadyABLInflowUFvPatchVectorField
    (
     const fvPatch&,
     const DimensionedField<vector, volMesh>&,
     const dictionary&
     );

    //- Construct by mapping given
    // unsteadyABLInflowUFvPatchVectorField onto a new patch
    unsteadyABLInflowUFvPatchVectorField
    (
     const unsteadyABLInflowUFvPatchVectorField&,
     const fvPatch&,
     const DimensionedField<vector, volMesh>&,
     const fvPatchFieldMapper&
     );

    //- Construct and return a clone
    virtual tmp<fvPatchVectorField> clone() const
    {
      return tmp<fvPatchVectorField>
	(
	 new unsteadyABLInflowUFvPatchVectorField(*this)
	 );
    }

    //- Construct as copy setting internal field reference
    unsteadyABLInflowUFvPatchVectorField
    (
     const unsteadyABLInflowUFvPatchVectorField&,
     const DimensionedField<vector, volMesh>&
     );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField> clone
    (
     const DimensionedField<vector, volMesh>& iF
     ) const
    {
      return tmp<fvPatchVectorField>
	(
	 new unsteadyABLInflowUFvPatchVectorField(*this, iF)
	 );
    }


    // Member functions
    // BC for vectorField is tricky. So using the function to return the vector field.
    tmp<vectorField> neutral(scalar M,vector flowDir,scalar MOL,scalar PBLH, scalar MBL) const;
    tmp<vectorField> veryStable(scalar M,vector flowDir,scalar MOL,scalar PBLH, scalar MBL) const;
    tmp<vectorField> stable(scalar M,vector flowDir,scalar MOL,scalar PBLH, scalar MBL) const;
    tmp<vectorField> unstable(scalar M,vector flowDir,scalar MOL,scalar PBLH, scalar MBL) const;
    //- Update the coefficients associated with the patch field
    using fixedValueFvPatchVectorField::updateCoeffs;
    virtual void updateCoeffs();
    //- Write
    virtual void write(Ostream&) const;
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
