const IOdictionary radProperties
(
 IOobject
 (
  "radiationProperties",    // dictionary name
  runTime.constant(),     // dict is found in "constant"
  mesh,                   // registry for the dict
  IOobject::MUST_READ,    // must exist, otherwise failure
  IOobject::NO_WRITE      // dict is only read by the solver                                                                 
  )
 );
const dictionary& solarDict = radProperties.subDict("solarLoadCoeffs");
List <List<scalar> > radProfileTable(solarDict.lookup("profileTable"));
scalarField radProfile(profileTable.size(),0.0);
scalarField radTimeProfile_(radProfileTable.size(),0.0);
scalarField radDirectswProfile_(radProfileTable.size(),0.0);
scalarField radDiffuseswProfile_(radProfileTable.size(),0.0);
forAll(radProfile,i)
{
  radTimeProfile_[i] = radProfileTable[i][0];
  radDirectswProfile_[i] = radProfileTable[i][1];
  radDiffuseswProfile_[i] = radProfileTable[i][2];
  Info<<" Radiation:"<<radTimeProfile_[i]<<"   "<<radDiffuseswProfile_[i]<<"   "<<endl;
}

