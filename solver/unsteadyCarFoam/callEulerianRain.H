#include "initializeRainFields.H"
#include "createTDFields.H"
Info << "Getting properties of air and water " << endl;
GET_parameters(temp,rhoa,mua,rhop);
Info << "Getting the rain phases and computing relevant quantities " << endl;
GET_fh(phases,fh,Rh,Vt);
for (label i=0;i<fh.size();i++)
  Info << "fh[" << i << "] = " << fh[i] << endl;	
fileName interFile(runTime.path()/"additionalData");
OFstream os(interFile);
for (label i=0;i<fh.size();i++)
  os << "fh" << (i+1) << " " << fh[i] << ";" << endl;
for (label i=0;i<Vt.size();i++)
  os << "Vt" << (i+1) << " " << Vt[i] << ";" << endl;	
os << "Rh" << " " << Rh.value() << ";" << endl;		
#include "createRainFields.H"
Info <<nl<< "Temperature: " << temp.value() << " K" << endl;
Info << "Air density: " << rhoa.value() << " kg/m3" << endl;
Info << "Air dynamic viscosity: " << mua.value() << " kg/m-s" << endl;
Info << "Water density: " << rhop.value() << " kg/m3" << endl;

