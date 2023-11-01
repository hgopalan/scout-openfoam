Foam::List<Foam::Pair<Foam::scalar> > CdMReM; // first is Cd and second is Re
Foam::List<Foam::Tuple2<Foam::scalar,Foam::Pair<Foam::scalar> > > allData; // fh, dia, vt

void GET_parameters (dimensionedScalar temp_, dimensionedScalar &rhoa_, dimensionedScalar &mua_, dimensionedScalar &rhop_)
{
	label i;
	
	scalar T_1_M[]={-18.0, -6.7, 0.0, 4.4, 15.6, 26.7, 37.8};
	scalar RHO_A_M[] = {1.38, 1.32, 1.293, 1.27, 1.22, 1.18, 1.13};
	scalar MU_A_M[] = {0.0157E-3, 0.0168E-3, 0.0171E-3, 0.0173E-3, 0.0179E-3, 0.0184E-3, 0.0190E-3};
	
	scalar T_2_M[] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0,
				20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0};
	scalar RHO_P_M[] = {999.87, 999.97, 1000.00, 999.97, 999.88, 999.73, 999.52, 999.27, 998.97, 998.62,
				998.23, 997.80, 997.33, 996.81, 996.26, 995.68, 994.00, 992.00};
				
	if (temp_.value()-scalar(273.15) < scalar(-18))
	{
		Info << "Air temperature is too low!";
	}
	else if (temp_.value()-scalar(273.15) >= scalar(37.8))
	{
		Info << "Air temperature is too high!";
	}
	else
	{
		for (i=0; i<=5; ++i)
		{
			if ( (T_1_M[i] <= temp_.value()-scalar(273.15)) && (temp_.value()-scalar(273.15) < T_1_M[i+1]) )
			{
				mua_.value() = MU_A_M[i] + ( ((MU_A_M[i+1] - MU_A_M[i])/(T_1_M[i+1] - T_1_M[i])) * (temp_.value()-scalar(273.15) - T_1_M[i]));
				rhoa_.value() = RHO_A_M[i] + ( ((RHO_A_M[i+1] - RHO_A_M[i])/(T_1_M[i+1] - T_1_M[i])) * (temp_.value()-scalar(273.15) - T_1_M[i]));
				break;
			}
		}
	}

	if (temp_.value()-scalar(273.15) <= scalar(0))
	{
		Info << "Air temperature is too low!";
	}
	else if (temp_.value()-scalar(273.15) >= scalar(40))
	{
		Info << "Air temperature is too high!";
	}
	else
	{
		for (i=0; i<=16; ++i)
		{
			if ( (T_2_M[i] <= temp_.value()-scalar(273.15)) && (temp_.value()-scalar(273.15) < T_2_M[i+1]) )
			{
				rhop_.value() = RHO_P_M[i] + ( ((RHO_P_M[i+1] - RHO_P_M[i])/(T_2_M[i+1] - T_2_M[i])) * (temp_.value()-scalar(273.15) - T_2_M[i]));
				break;
			}
		}
	}		
}

volScalarField GET_CdRe (volScalarField Reynolds)
{
	volScalarField CdRe = Reynolds; //just initialization
	forAll (Reynolds, celli) 
	{
		label i;
	
		
		
		if (Reynolds[celli] < scalar(1.8))
		{
			i = 0;
			CdRe[celli] = ( CdMReM[i].first() + (((CdMReM[i+1].first() - CdMReM[i].first())/(CdMReM[i+1].second() - CdMReM[i].second()))*(Reynolds[celli] - CdMReM[i].second())) ) * Reynolds[celli];
		}
		else if (Reynolds[celli] >= scalar(3549))
		{
			i = 32;
			CdRe[celli] = ( CdMReM[i].first() + (((CdMReM[i+1].first() - CdMReM[i].first())/(CdMReM[i+1].second() - CdMReM[i].second()))*(Reynolds[celli] - CdMReM[i].second())) ) * Reynolds[celli];
		}
		else
		{
			for (i=0; i<=32; ++i)
			{
				if ( (CdMReM[i].second() <= Reynolds[celli]) && (Reynolds[celli] < CdMReM[i+1].second()) )
				{
					CdRe[celli] = ( CdMReM[i].first() + (((CdMReM[i+1].first() - CdMReM[i].first())/(CdMReM[i+1].second() - CdMReM[i].second()))*(Reynolds[celli] - CdMReM[i].second())) ) * Reynolds[celli];
					break;
				}
			}
		}

	}
	
	return (CdRe);
}

/*
scalar best_result (scalar d,scalar rh)
{
	scalar a,A,n,p,res,h;
	A = 1.3;
	n = 2.25;
	p = 0.232;
	
	a = A*Foam::pow(rh,p);
	d*=1000;
	h = d/a;
	res = 1 - Foam::exp(-Foam::pow((h),n));
	return res;	
}
*/

scalar rainDistribution (scalar d, scalar rh, word dist="Best")
{
	scalar res=0.0;
	
	if(dist == "Best")
	{
		scalar a,A,n,p,h;
		A = 1.3;
		n = 2.25;
		p = 0.232;

		a = A*Foam::pow(rh,p);
		d*=1000;
		h = d/a;
		res = 1 - Foam::exp(-Foam::pow(h,n));
	}
	return res;
}

void GET_fh (List<scalar> dia_, List<scalar> &fh_, dimensionedScalar Rh_, List<scalar> &Vt_)
{
	scalar dp, dn, fhSum = 0;
	label i,j;
	
	for(i=1;i<allData.size()-1;i++)
	{
		dp = allData[i].second().first() + (allData[i+1].second().first()-allData[i].second().first())/2;
		dn = allData[i].second().first() - (allData[i].second().first()-allData[i-1].second().first())/2;
		allData[i].first() = rainDistribution(dp,Rh_.value()) - rainDistribution(dn,Rh_.value());
		fhSum+=allData[i].first()*allData[i].second().second();
	}
	
	dp = allData[0].second().first() + (allData[1].second().first()-allData[0].second().first())/2;
	dn = allData[allData.size()-1].second().first() - (allData[allData.size()-1].second().first()-allData[allData.size()-2].second().first())/2;

	allData[0].first() = rainDistribution(dp,Rh_.value());
	allData[allData.size()-1].first() = 1-rainDistribution(dn,Rh_.value());
	
	for(i=0;i<dia_.size();i++)
	{
		for (j=0;j<allData.size();j++)
			if(dia_[i]==allData[j].second().first())
			{
			  fh_[i]=(allData[j].second().second()*allData[j].first());///fhSum;
				Info<<"Harish :"<<fh_[i]<<"  "<<fhSum<<endl;
				Vt_[i]=allData[j].second().second();
				break;			
			}
	}
}
