{
  // U
  List <List<scalar> > cCRSProfileTableU(ABLDict.lookup("cCRSProfileTableU"));
  scalar cCRSStartTime (ABLDict.lookupOrDefault<scalar>("cCRSStartTime",13));
  scalar timeValue=U.time().value()+cCRSStartTime*3600;
  scalar cCRScounter=int(timeValue/3600);
  scalar fraction=timeValue/3600-int(timeValue/3600);
  scalar cCRSProfilesize (ABLDict.lookupOrDefault<scalar>("cCRSProfileSize",100));
  //Info<<"CCRS Profile Size:"<<cCRSProfilesize<<endl;
  scalarField cCRSdummyscalar(cCRSProfilesize,0.0);
  scalarIOList cCRStimeprofile_
    (
     IOobject
     (
      "cCRStimeProfile",
      runTime.constant(),
      mesh,
      IOobject::NO_READ,    // must exist, otherwise failure
      IOobject::NO_WRITE      // dict is only read by the solver
      ),
     cCRSdummyscalar
     );
  scalarIOList cCRSlowerprofile_
   (
     IOobject
     (
      "cCRSlowerProfile",
      runTime.constant(),
      mesh,
      IOobject::NO_READ,    // must exist, otherwise failure
      IOobject::NO_WRITE      // dict is only read by the solver
      ),
     cCRSdummyscalar
    );
 scalarIOList cCRSupperprofile_
   (
     IOobject
     (
      "cCRSupperProfile",
      runTime.constant(),
      mesh,
      IOobject::NO_READ,    // must exist, otherwise failure
      IOobject::NO_WRITE      // dict is only read by the solver
      ),
     cCRSdummyscalar
    );
 scalarField ccrstimevalues(cCRStimeprofile_.size(),1e20);
 scalarField ccrslvalues(cCRStimeprofile_.size(),0.0);
 scalarField ccrsuvalues(cCRStimeprofile_.size(),0.0);
 //Info<<"CCRSSSSS:"<<cCRScounter<<endl;
 //Info<<"Size:"<<cCRStimeprofile_.size()<<endl;
 forAll(cCRStimeprofile_,i)
   {
     cCRStimeprofile_[i]=cCRSProfileTableU[i][0];
     cCRSlowerprofile_[i]=cCRSProfileTableU[i][cCRScounter+1];
     cCRSupperprofile_[i]=cCRSProfileTableU[i][cCRScounter+2];
     //Info<<i<<"  "<<cCRStimeprofile_[i]<<"   "<<cCRSlowerprofile_[i]<<"  "<<cCRSupperprofile_[i]<<endl;
   }
 //Info<<"U Wind:"<<endl;
 forAll(ccrstimevalues,i)
   {
     ccrstimevalues[i] = cCRStimeprofile_[i];
     ccrslvalues[i] = cCRSlowerprofile_[i];
     ccrsuvalues[i] = cCRSupperprofile_[i];
     //Info<<ccrstimevalues[i]<<"   "<<ccrslvalues[i]<<"  "<<ccrsuvalues[i]<<endl;
   }
 forAll(zProfile_,i)
   {
     scalar z=zProfile_[i];
     scalar lowerValue=interpolateXY(z,ccrstimevalues,ccrslvalues); 
     scalar upperValue=interpolateXY(z,ccrstimevalues,ccrsuvalues);
     if(mag(lowerValue)>0 && mag(upperValue)<1e-5)
       {
	vertU[i]=lowerValue;
	fraction=0;
       }
       else if(z>ccrstimevalues[0])
       vertU[i]=(1-fraction)*lowerValue+fraction*upperValue;
     //Info<<zProfile_[i]<<"  "<<vertU[i]<<endl;
   }
  //V
    List <List<scalar> > cCRSProfileTableV(ABLDict.lookup("cCRSProfileTableV"));
forAll(cCRStimeprofile_,i)
   {
     cCRStimeprofile_[i]=cCRSProfileTableV[i][0];
     cCRSlowerprofile_[i]=cCRSProfileTableV[i][cCRScounter+1];
     cCRSupperprofile_[i]=cCRSProfileTableV[i][cCRScounter+2];
   }
 forAll(ccrstimevalues,i)
   {
     ccrstimevalues[i] = cCRStimeprofile_[i];
     ccrslvalues[i] = cCRSlowerprofile_[i];
     ccrsuvalues[i] = cCRSupperprofile_[i];
     //Info<<ccrstimevalues[i]<<"   "<<ccrslvalues[i]<<"  "<<ccrsuvalues[i]<<endl;
   }
 forAll(zProfile_,i)
   {
     scalar z=zProfile_[i];
     scalar lowerValue=interpolateXY(z,ccrstimevalues,ccrslvalues); 
     scalar upperValue=interpolateXY(z,ccrstimevalues,ccrsuvalues);
     if(mag(lowerValue)>0 && mag(upperValue)<1e-5)
       {
        vertV[i]=lowerValue;
       }
       else if(z>ccrstimevalues[0])
	 vertV[i]=(1-fraction)*lowerValue+fraction*upperValue; 
     //Info<<zProfile_[i]<<"  "<<vertU[i]<<endl;
   }
  //theta
    List <List<scalar> > cCRSProfileTableT(ABLDict.lookup("cCRSProfileTableT"));
forAll(cCRStimeprofile_,i)
   {
     cCRStimeprofile_[i]=cCRSProfileTableT[i][0];
     cCRSlowerprofile_[i]=cCRSProfileTableT[i][cCRScounter+1];
     cCRSupperprofile_[i]=cCRSProfileTableT[i][cCRScounter+2];
   }
 forAll(ccrstimevalues,i)
   {
     ccrstimevalues[i] = cCRStimeprofile_[i];
     ccrslvalues[i] = cCRSlowerprofile_[i];
     ccrsuvalues[i] = cCRSupperprofile_[i];
     //Info<<ccrstimevalues[i]<<"   "<<ccrslvalues[i]<<"  "<<ccrsuvalues[i]<<endl;
   }
 forAll(zProfile_,i)
   {
     scalar z=zProfile_[i];
     scalar lowerValue=interpolateXY(z,ccrstimevalues,ccrslvalues); 
     scalar upperValue=interpolateXY(z,ccrstimevalues,ccrsuvalues);
     if(mag(lowerValue)>0 && mag(upperValue)<1e-5)
       {
        vertT[i]=lowerValue;
       }
       else if(z>ccrstimevalues[0])
	 vertT[i]=(1-fraction)*lowerValue+fraction*upperValue; 
     //Info<<zProfile_[i]<<"  "<<vertU[i]<<endl;
   }
  //qv
    List <List<scalar> > cCRSProfileTableq(ABLDict.lookup("cCRSProfileTableq"));
forAll(cCRStimeprofile_,i)
   {
     cCRStimeprofile_[i]=cCRSProfileTableq[i][0];
     cCRSlowerprofile_[i]=cCRSProfileTableq[i][cCRScounter+1];
     cCRSupperprofile_[i]=cCRSProfileTableq[i][cCRScounter+2];
   }
 forAll(ccrstimevalues,i)
   {
     ccrstimevalues[i] = cCRStimeprofile_[i];
     ccrslvalues[i] = cCRSlowerprofile_[i];
     ccrsuvalues[i] = cCRSupperprofile_[i];
     //Info<<ccrstimevalues[i]<<"   "<<ccrslvalues[i]<<"  "<<ccrsuvalues[i]<<endl;
   }
 forAll(zProfile_,i)
   {
     scalar z=zProfile_[i];
     scalar lowerValue=interpolateXY(z,ccrstimevalues,ccrslvalues); 
     scalar upperValue=interpolateXY(z,ccrstimevalues,ccrsuvalues);
     if(mag(lowerValue)>0 && mag(upperValue)<1e-5)
       {
        vertqV[i]=lowerValue;
       }
       else if(z>ccrstimevalues[0])
	 vertqV[i]=(1-fraction)*lowerValue+fraction*upperValue;
     /* if(zLoc[i]<1.2*PBLH)
     Info<<zProfile_[i]<<"  "<<vertU[i]<<" "<<vertV[i]<<" "
     <<vertT[i]<<" "<<vertqV[i]<<endl; */
   }
}
