PtrList<volScalarField> ps;
ps.setSize(scalarTable.size());
std::stringstream converter; 	
for (int phase_no = 0; phase_no < scalarTable.size(); phase_no++)
  {
    string phase_Uname = "ps";
    string phase_no_str;
    std::stringstream convert;
    convert << phase_no+1;
    phase_no_str = convert.str();
    phase_Uname.append(phase_no_str);
    ps.set(phase_no,new volScalarField
	   (
	    IOobject
	    (
	     phase_Uname,
	     runTime.timeName(),
	     mesh,
	     IOobject::READ_IF_PRESENT,
	     IOobject::AUTO_WRITE
	     ),
	    0+0.0*mag(U),
	    zeroGradientFvPatchField<scalar>::typeName
	    ));
  }

	

for (int i=0; i< scalarTable.size();++i)                                                                              
  ps[i].write();   



