/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copypdright (C) 2011-2015 OpenFOAM Foundation
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

#include "unsteadyABLInflowUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interpolateXY.H"
#include "scalarIOList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  unsteadyABLInflowUFvPatchVectorField::
  unsteadyABLInflowUFvPatchVectorField
  (
   const fvPatch& p,
   const DimensionedField<vector, volMesh>& iF
   )
    :
    fixedValueFvPatchVectorField(p, iF),
    kappa_(0.41),
    ZRef_(10.0),
    z0_(0.1),
    zd_(0.0)

  {}


  unsteadyABLInflowUFvPatchVectorField::
  unsteadyABLInflowUFvPatchVectorField
  (
   const fvPatch& p,
   const DimensionedField<vector, volMesh>& iF,
   const dictionary& dict
   )
    :
    fixedValueFvPatchVectorField(p, iF),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    ZRef_(dict.lookupOrDefault<scalar>("ZRef",10.0)),   
    z0_(dict.lookupOrDefault<scalar>("z0", 0.1)),
    zd_(dict.lookupOrDefault<scalar>("zd", 0.0))    
  {
    if (dict.found("value"))
      {
	fvPatchField<vector>::operator=
	  (
	   vectorField("value", dict, p.size())
	   );
      }
    else
      {
	updateCoeffs();
      } 
  }


  unsteadyABLInflowUFvPatchVectorField::
  unsteadyABLInflowUFvPatchVectorField
  (
   const unsteadyABLInflowUFvPatchVectorField& pvf,
   const fvPatch& p,
   const DimensionedField<vector, volMesh>& iF,
   const fvPatchFieldMapper& mapper
   )
    :
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    kappa_(pvf.kappa_),
    ZRef_(pvf.ZRef_),    
    z0_(pvf.z0_),
    zd_(pvf.zd_)
  {}


  unsteadyABLInflowUFvPatchVectorField::
  unsteadyABLInflowUFvPatchVectorField
  (
   const unsteadyABLInflowUFvPatchVectorField& pvf,
   const DimensionedField<vector, volMesh>& iF
   )
    :
    fixedValueFvPatchVectorField(pvf, iF),
    kappa_(pvf.kappa_),
    ZRef_(pvf.ZRef_),        
    z0_(pvf.z0_),
    zd_(pvf.zd_)    
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


  void unsteadyABLInflowUFvPatchVectorField::updateCoeffs()
  {
    if (this->updated())
      {
	return;
      }
    const volScalarField p_rgh_=this->db().lookupObject<volScalarField>("p_rgh"); 
    // The reference wind speed is read from constant/ABLDict. 
    // Default: (Time U V T MOL)
    const scalarIOList&  Tvalues=this->db().lookupObject<scalarIOList>("timeProfile");
    const scalarIOList& Uvalues=this->db().lookupObject<scalarIOList>("uProfile");
    const scalarIOList& Vvalues=this->db().lookupObject<scalarIOList>("vProfile");
    const scalarIOList& MOvalues=this->db().lookupObject<scalarIOList>("MOProfile");
    const scalarIOList& PBLHvalues=this->db().lookupObject<scalarIOList>("PBLHProfile");
    const scalarIOList& MBLvalues=this->db().lookupObject<scalarIOList>("MBLProfile");
    scalarField timevalues(Tvalues.size(),0.0);
    scalarField uvalues(Tvalues.size(),0.0);
    scalarField vvalues(Tvalues.size(),0.0);
    scalarField movalues(Tvalues.size(),0.0);
    scalarField pblhvalues(Tvalues.size(),0.0);
    scalarField mblvalues(Tvalues.size(),0.0);
    forAll(Tvalues,i)
      {
	timevalues[i] = Tvalues[i];
	uvalues[i] = Uvalues[i];
	vvalues[i] = Vvalues[i];
	movalues[i] = MOvalues[i];
	pblhvalues[i] = PBLHvalues[i];
	mblvalues[i] = MBLvalues[i];
      }
    scalar u1=interpolateXY(p_rgh_.time().value(),timevalues,uvalues);
    scalar v1=interpolateXY(p_rgh_.time().value(),timevalues,vvalues);
    scalar M=sqrt(sqr(u1)+sqr(v1));
    scalar MOL=interpolateXY(p_rgh_.time().value(),timevalues,movalues);
    if(mag(MOL)>25)
      MOL=MOL;
    else if(neg(MOL) && MOL> -25)
      MOL=-25;
    else if(pos(MOL) && mag(MOL)<25)
      MOL=25;   
    scalar PBLH=interpolateXY(p_rgh_.time().value(),timevalues,pblhvalues);
    scalar MBL=interpolateXY(p_rgh_.time().value(),timevalues,mblvalues);
    const vector flowDir(u1/M,v1/M,0);
    // Neutral
    if(mag(MOL)>500)
      {
	vectorField::operator=(neutral(M,flowDir,MOL,PBLH,MBL)); 
	fixedValueFvPatchVectorField::updateCoeffs();
      }
     else if(pos(MOL) && MOL<2.5)
     {
	vectorField::operator=(veryStable(M,flowDir,MOL,PBLH,MBL)); 
	fixedValueFvPatchVectorField::updateCoeffs();
	 }
	 else if(pos(MOL))
    {
	vectorField::operator=(stable(M,flowDir,MOL,PBLH,MBL)); 
	fixedValueFvPatchVectorField::updateCoeffs();
	 }
	 else 
    {
	vectorField::operator=(unstable(M,flowDir,MOL,PBLH,MBL)); 
	fixedValueFvPatchVectorField::updateCoeffs();
	 }
  }
  tmp<vectorField> unsteadyABLInflowUFvPatchVectorField::neutral(scalar M,vector flowDir,scalar MOL, scalar PBLH, scalar MBL) const
  {
    // If not inflow just set to zeroGradient
    scalar fluxdir=gSum(patch().nf() & flowDir);
    if(fluxdir>0 || patch().name()=="upper")
      return this->patchInternalField();
    // Compute u*
    scalar denominator=log(ZRef_/z0_)+ZRef_/MBL-ZRef_/PBLH*0.5*ZRef_/MBL;
    scalar ustar=M*kappa_/denominator;
    // Compute z and restrict height to PBLH 
    scalarField z_=patch().Cf().component(vector::Z)-zd_+z0_;
    z_=min(z_,0.9*PBLH);
    // Compute WindSpeed and return the value
    return (flowDir*ustar/kappa_*(log(z_/z0_)+z_/MBL-z_/PBLH*0.5*z_/MBL));
  }
 tmp<vectorField> unsteadyABLInflowUFvPatchVectorField::veryStable(scalar M,vector flowDir,scalar MOL, scalar PBLH, scalar MBL) const
  {
    // If not inflow just set to zeroGradient
    scalar fluxdir=gSum(patch().nf() & flowDir);
    if(fluxdir>0 || patch().name()=="upper")
      return this->patchInternalField();
    // Compute u*
    scalar denominator=6*log(ZRef_/z0_)-5*ZRef_/PBLH+ZRef_/MBL-ZRef_/PBLH*0.5*ZRef_/MBL;
    scalar ustar=M*kappa_/denominator;
    // Compute z and restrict height to PBLH 
    scalarField z_=patch().Cf().component(vector::Z)-zd_+z0_;
    z_=min(z_,0.9*PBLH);
    // Compute WindSpeed and return the value
    return (flowDir*ustar/kappa_*(6*log(z_/z0_)-5*z_/PBLH+z_/MBL-z_/PBLH*0.5*z_/MBL));
  }
tmp<vectorField> unsteadyABLInflowUFvPatchVectorField::stable(scalar M,vector flowDir,scalar MOL, scalar PBLH, scalar MBL) const
  {
    // If not inflow just set to zeroGradient
    scalar fluxdir=gSum(patch().nf() & flowDir);
    if(fluxdir>0 || patch().name()=="upper")
      return this->patchInternalField();
    // Compute u*
    scalar denominator=log(ZRef_/z0_)+5*ZRef_/MOL*(1-0.5*ZRef_/PBLH)+ZRef_/MBL-ZRef_/PBLH*0.5*ZRef_/MBL;
    scalar ustar=M*kappa_/denominator;
    // Compute z and restrict height to PBLH 
    scalarField z_=patch().Cf().component(vector::Z)-zd_+z0_;
    z_=min(z_,0.9*PBLH);
    // Compute WindSpeed and return the value
    return (flowDir*ustar/kappa_*(log(z_/z0_)+5*z_/MOL*(1-0.5*z_/PBLH)+z_/MBL-z_/PBLH*0.5*z_/MBL));
  }
  
tmp<vectorField> unsteadyABLInflowUFvPatchVectorField::unstable(scalar M,vector flowDir,scalar MOL, scalar PBLH, scalar MBL) const
  {
    // If not inflow just set to zeroGradient
    scalar fluxdir=gSum(patch().nf() & flowDir);
    if(patch().name()=="upper")
      {
	return this->patchInternalField();
      }
    // Compute u*
    scalar x=pow(1-12*ZRef_/MOL,1.0/3.0);
    scalar psiM=1.5*log((1+x+sqr(x))/3.0)+3.14/sqrt(3.0)-sqrt(3.0)*atan((1+2*x)/sqrt(3.0));
    scalar denominator=log(ZRef_/z0_)-psiM+ZRef_/PBLH*(1+(pow(1-12*ZRef_/MOL,1.0/3.0)-1)/(8*ZRef_/MOL));
    denominator=denominator+ZRef_/MBL-ZRef_/PBLH*0.5*ZRef_/MBL;
    scalar ustar=M*kappa_/denominator;
    // On water friction velocity is a mess ; so limting it 
    // Compute z and restrict height to PBLH 
    scalarField z_=patch().Cf().component(vector::Z)-zd_+z0_;
    z_=min(z_,0.9*PBLH);
    // Compute WindSpeed and return the value
    scalarField xF=pow(1-12*z_/MOL,1.0/3.0);
    scalarField psiMF=1.5*log((1+xF+sqr(xF))/3.0)+3.14/sqrt(3.0)-sqrt(3.0)*atan((1+2*xF)/sqrt(3.0));
    scalarField windSpeed=log(z_/z0_)-psiMF+z_/PBLH*(1+(pow(1-12*z_/MOL,1.0/3.0)-1)/(8*z_/MOL));
    windSpeed=windSpeed+z_/MBL-z_/PBLH*0.5*z_/MBL;
    if(fluxdir>0)
      {
	//	scalarField numericalFlux= (this->patchInternalField() & patch().Sf());
	//	scalarField analyticalFlux=((flowDir*ustar/kappa_*windSpeed) & patch().Sf());
	//	vectorField acceleration=(numericalFlux - analyticalFlux)*flowDir/patch().magSf();
	//	return ((flowDir*ustar/kappa_*windSpeed)+acceleration);
	return this->patchInternalField();
      }
    return (flowDir*ustar/kappa_*windSpeed);
  }
  void unsteadyABLInflowUFvPatchVectorField::write(Ostream& os) const
  {
    fvPatchVectorField::write(os);
    os.writeKeyword("kappa")<< kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("ZRef")<< ZRef_ << token::END_STATEMENT << nl;            
    os.writeKeyword("z0")<< z0_ << token::END_STATEMENT << nl;
    os.writeKeyword("zd")<< zd_ << token::END_STATEMENT << nl;    
    writeEntry("value", os);
  }


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  makePatchTypeField
  (
   fvPatchVectorField,
   unsteadyABLInflowUFvPatchVectorField
   );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam



// ************************************************************************* //
