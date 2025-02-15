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
  ! Wall boundary condition using Monin-Obukhov Similarity Theory              !
  !----------------------------------------------------------------------------!                                                         
  MOST requires that nut dU/dn = ustar^2. Fixing the velocity gradient at the wall
  is complicated in OpenFOAM as it is a tensor quantity. Instead we use an indirect
  implementation as follows. 
  
  Step 1: Calculat the cell-center velocity parallel to the surface 
  \vec{U}=\vec{U_c}-\dot{U_c,nhatV}nhatV
  U_c is the cell-center velocity vector. nhatV is the normal vector between cell and
  face center. 
  Step 2: Calculate ustar from U using the following expression
  |U| = ustar/kappa * [ log ((nc+z0)/z0)-psiM((nc+z0)/L)], where nc is distance from
  face to cell center. 
  Step 3: Calculate nut = ustar * kappa * (nc+z0)/phiM((nc+z0)/L)
  Step 4: Replace dU/dn with (U_c - U_w)/nc. U_w is the value we fix at the wall face-
  center to ensure that nut dU/dn = ustar^2. It should be noted that U_w is not a 
  physical value but a numerical value to ensure consistency. 
  Step 5: Calculate dU = ustar*nc/(kappa*(nc+z0))*phiM
  Step 6: \vec{U_w} = \vec{U}*(1-dU/|U|)
  The above formulation is valid for both flat and curved surfaces. 

  References: 
  (HG2018) Gopalan, H., "Evaluation of Wray-Agarwal turbulence model for simulation of 
  neutral and non-neutral atmospheric boundary layers",JWEIA, 2018: 182, 322-329.

  Log File:
  1. 28/3/2019: Adopted for OF-v1812. 

  Authors:
  Harish Gopalan (gopalanh@ihpc.a-star.edu.sg)
*/ 

/* Usage 
   terrain
   {
   type ABLMOWallU;
   z0 uniform 1.000000; // Roughness length in m. Default is 0.1 m. 
   value uniform (0 0 0); // This is a placeholder. 
   }
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef ABLMOWallUFvPatchVectorField_H
#define ABLMOWallUFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  namespace incompressible
  {
    /*---------------------------------------------------------------------------*\
      Class ABLMOWallUFvPatchVectorField Declaration
      \*---------------------------------------------------------------------------*/

    class ABLMOWallUFvPatchVectorField
      :
      public fixedValueFvPatchVectorField

    {
      scalar kappa_; // von-Karman constant, 0.41 
      scalar z0_; // Roughness length in m. Should make it a scalar in future....

    public:

      //- Runtime type information
      TypeName("ABLMOWallU");


      // Constructors

      //- Construct from patch and internal field
      ABLMOWallUFvPatchVectorField
      (
       const fvPatch&,
       const DimensionedField<vector, volMesh>&
       );

      //- Construct from patch, internal field and dictionary
      ABLMOWallUFvPatchVectorField
      (
       const fvPatch&,
       const DimensionedField<vector, volMesh>&,
       const dictionary&
       );

      //- Construct by mapping given
      // ABLMOWallUFvPatchVectorField onto a new patch
      ABLMOWallUFvPatchVectorField
      (
       const ABLMOWallUFvPatchVectorField&,
       const fvPatch&,
       const DimensionedField<vector, volMesh>&,
       const fvPatchFieldMapper&
       );

      //- Construct and return a clone
      virtual tmp<fvPatchVectorField> clone() const
      {
	return tmp<fvPatchVectorField>
	  (
	   new ABLMOWallUFvPatchVectorField(*this)
	   );
      }

      //- Construct as copy setting internal field reference
      ABLMOWallUFvPatchVectorField
      (
       const ABLMOWallUFvPatchVectorField&,
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
	   new ABLMOWallUFvPatchVectorField(*this, iF)
	   );
      }


      // Member functions
      tmp<vectorField> neutral(scalar M,vector flowDir,scalar MOL,scalar PBLH, scalar MBL) const;
      tmp<vectorField> veryStable(scalar M,vector flowDir,scalar MOL,scalar PBLH, scalar MBL) const;
      tmp<vectorField> stable(scalar M,vector flowDir,scalar MOL,scalar PBLH, scalar MBL) const;
      tmp<vectorField> unstable(scalar M,vector flowDir,scalar MOL,scalar PBLH, scalar MBL) const;
      // Mapping functions

      //- Map (and resize as needed) from self given a mapping object
      virtual void autoMap
      (
       const fvPatchFieldMapper&
       );

      //- Reverse map the given fvPatchField onto this fvPatchField
      virtual void rmap
      (
       const fvPatchVectorField&,
       const labelList&
       );

      //- Update the coefficients associated with the patch field
      using fixedValueFvPatchVectorField::updateCoeffs;
      virtual void updateCoeffs();
      //- Write
      virtual void write(Ostream&) const;
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  }
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
