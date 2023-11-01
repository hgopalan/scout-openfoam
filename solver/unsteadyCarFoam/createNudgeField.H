
volVectorField U_avg=0*U;
//volVectorField Unudge=0*U;
//volScalarField nutnudge=0*nut;
scalar verticalprofilesize_=100;
scalarField verticaldummyscalar(verticalprofilesize_,0.0);
scalarIOList verticalzprofile_
(
 IOobject
 (
  "verticalzProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  ),
 verticaldummyscalar
 );
scalarIOList verticaluprofile_
(
 IOobject
 (
  "verticaluProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  ),
 verticaldummyscalar
 );
scalarIOList verticalvprofile_
(
 IOobject
 (
  "verticalvProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  ),
 verticaldummyscalar
 );
scalarIOList verticalTprofile_
(
 IOobject
 (
  "verticalTProfile",
  runTime.constant(),
  mesh,
  IOobject::NO_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver
  ),
 verticaldummyscalar
 );

scalarField verticallevels(verticalzprofile_.size(),0.0);
scalarField verticalzvalues(verticalzprofile_.size(),0.0);
scalarField verticalURefvalues(verticalzprofile_.size(),0.0);
scalarField verticalVRefvalues(verticalzprofile_.size(),0.0);
scalarField verticalTRefvalues(verticalzprofile_.size(),0.0);
scalarField verticalnutRefvalues(verticalzprofile_.size(),0.0);

// Experimental 
volScalarField cz=(mesh.C()).component(vector::Z);
scalar zmin=(min(cz)).value();
scalar zmax=(max(cz)).value();
reduce(zmin,minOp<scalar>());
reduce(zmax,maxOp<scalar>());
Info<<"Limits:"<<zmin<<"  "<<zmax<<endl;

forAll(verticalzprofile_,i)
{
  verticaldummyscalar[i]=i;
  verticalzprofile_[i]=zmin+(i+1)*(zmax-zmin)/verticalzprofile_.size();
}

forAll(verticalzprofile_,i)
{
  verticalzvalues[i] = verticalzprofile_[i];
  verticallevels[i]=verticaldummyscalar[i];
}

volScalarField verticallevelMarker=0*T;
forAll(verticallevelMarker,cellI)
{
  verticallevelMarker[cellI]=int(interpolateXY(cz[cellI],verticalzvalues,verticallevels));
}

scalarField uavg(verticalzprofile_.size(),0.0);
scalarField vavg(verticalzprofile_.size(),0.0);
scalarField Tavg(verticalzprofile_.size(),0.0);


