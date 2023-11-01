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
  ! Temperature Inflow Conditions for Atmospheric Boundary Layer               !
  !----------------------------------------------------------------------------!                                                         
  This module implements the inflow boundary conditions for turbulent viscosity
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
   type unsteadyABLInflowQv;
   ZRef 15.000000; // Default value is 2.0 m 
   z0 1.000000; // Default value is 0.1 m 
   value uniform 300.15; // This is a placeholder. Use values in Kelvin. 
   }
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef unsteadyABLInflowQvFvPatchScalarField_H
#define unsteadyABLInflowQvFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
      Class unsteadyABLInflowQvFvPatchScalarField Declaration
      \*---------------------------------------------------------------------------*/

    class unsteadyABLInflowQvFvPatchScalarField
      :
      public fixedValueFvPatchScalarField

    {

    public:

      //- Runtime type information
      TypeName("unsteadyABLInflowQv");


      // Constructors

      //- Construct from patch and internal field
      unsteadyABLInflowQvFvPatchScalarField
      (
       const fvPatch&,
       const DimensionedField<scalar, volMesh>&
       );

      //- Construct from patch, internal field and dictionary
      unsteadyABLInflowQvFvPatchScalarField
      (
       const fvPatch&,
       const DimensionedField<scalar, volMesh>&,
       const dictionary&
       );

      //- Construct by mapping given
      // unsteadyABLInflowQvFvPatchScalarField onto a new patch
      unsteadyABLInflowQvFvPatchScalarField
      (
       const unsteadyABLInflowQvFvPatchScalarField&,
       const fvPatch&,
       const DimensionedField<scalar, volMesh>&,
       const fvPatchFieldMapper&
       );
      //- Construct as copy
      unsteadyABLInflowQvFvPatchScalarField
      (
       const unsteadyABLInflowQvFvPatchScalarField&
       );
      //- Construct and return a clone
      virtual tmp<fvPatchScalarField> clone() const
      {
	return tmp<fvPatchScalarField>
	  (
	   new unsteadyABLInflowQvFvPatchScalarField(*this)
	   );
      }

      //- Construct as copy setting internal field reference
      unsteadyABLInflowQvFvPatchScalarField
      (
       const unsteadyABLInflowQvFvPatchScalarField&,
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
	   new unsteadyABLInflowQvFvPatchScalarField(*this, iF)
	   );
      }


      // Member functions

      // Return access to z0
      // Returns temperature profile 
      tmp<scalarField> TRef(const vector& flowdir) const;     
      //- Update the coefficients associated with the patch field
      virtual void updateCoeffs();
      //- Write
      virtual void write(Ostream&) const;
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
