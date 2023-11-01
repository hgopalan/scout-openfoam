volScalarField::Boundary& pBf = p_rgh.boundaryFieldRef();
const volVectorField::Boundary& UBf = U.boundaryField();
const surfaceScalarField::Boundary& phiHbyABf = phiHbyA.boundaryField();
const surfaceScalarField::Boundary& rhorAUBf =rAUf.boundaryField();
const surfaceVectorField::Boundary& SfBf =mesh.Sf().boundaryField();
const surfaceScalarField::Boundary& magSfBf=mesh.magSf().boundaryField();
forAll(pBf, patchi)
{
  if (isA<fixedFluxPressureFvPatchScalarField>(pBf[patchi]))
    {
      refCast<fixedFluxPressureFvPatchScalarField>
	(
	 pBf[patchi]
	 ).updateSnGrad
	(
	 (
	  phiHbyABf[patchi]
	  - rhok.boundaryField()[patchi]
	  *MRF.relative(SfBf[patchi] & UBf[patchi], patchi)
	  )
	 /(magSfBf[patchi]*rhorAUBf[patchi])
	 );
    }
}
