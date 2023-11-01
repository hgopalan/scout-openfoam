/*
//--------------------------------------------------------------------------------------
		// Matrices containing Reynolds numbers and the corresponding drag coefficient values
		// according to Gunn & Kinzer
		Foam::scalar Re_M[]={1.80, 9.61, 23.4, 43.2, 68.7, 98.9, 134.0, 175.0, 220.0, 269.0,
			  372.0, 483.0, 603.0, 731.0, 866.0, 1013.0, 1164.0, 1313.0, 1461.0, 1613.0,
			  1764.0, 1915.0, 2066.0, 2211.0, 2357.0, 2500.0, 2636.0, 2772.0, 2905.0, 3033.0,
			  3164.0, 3293.0, 3423.0, 3549.0}; 
		Foam::scalar Cd_M[]={15.0, 4.2, 2.4, 1.66, 1.28, 1.07, 0.926, 0.815, 0.729, 0.671,
			  0.607, 0.570, 0.545, 0.528, 0.517, 0.504, 0.495, 0.494, 0.498, 0.503,
			  0.511, 0.520, 0.529, 0.544, 0.559, 0.575, 0.594, 0.615, 0.635, 0.660,
			  0.681, 0.700, 0.727, 0.751};
		//--------------------------------------------------------------------------------------

	Foam::scalar allDia[]={0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.0012,0.0014,0.0016,
			0.0018,0.002,0.003,0.004,0.005,0.006};
	Foam::scalar allFh[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Foam::scalar allVt[]={1.17,1.62,2.06,2.47,2.87,3.27,3.67,4.03,4.64,5.17,5.65,6.09,6.49,8.06,8.83,9.09,9.18};
*/

//using namespace Foam;

//Foam::List<Foam::Pair<Foam::scalar> > CdMReM; // first is Cd and second is Re
//CdMReM<Foam::Pair<Foam::scalar> >.setSize(34);
CdMReM.setSize(34);

//CdMReM<Foam::Pair<Foam::scalar> >[0] = Foam::Pair<Foam::scalar>(15,1.8);
CdMReM[0] = Foam::Pair<Foam::scalar>(15,1.8);
CdMReM[1] = Foam::Pair<Foam::scalar>(4.2,9.61);
CdMReM[2] = Foam::Pair<Foam::scalar>(2.4,23.4);
CdMReM[3] = Foam::Pair<Foam::scalar>(1.66,43.2);
CdMReM[4] = Foam::Pair<Foam::scalar>(1.28,68.7);
CdMReM[5] = Foam::Pair<Foam::scalar>(1.07,98.9);
CdMReM[6] = Foam::Pair<Foam::scalar>(0.926,134);
CdMReM[7] = Foam::Pair<Foam::scalar>(0.815,175);
CdMReM[8] = Foam::Pair<Foam::scalar>(0.729,220);
CdMReM[9] = Foam::Pair<Foam::scalar>(0.671,269);
CdMReM[10] = Foam::Pair<Foam::scalar>(0.607,372);
CdMReM[11] = Foam::Pair<Foam::scalar>(0.57,483);
CdMReM[12] = Foam::Pair<Foam::scalar>(0.545,603);
CdMReM[13] = Foam::Pair<Foam::scalar>(0.528,731);
CdMReM[14] = Foam::Pair<Foam::scalar>(0.517,866);
CdMReM[15] = Foam::Pair<Foam::scalar>(0.504,1013);
CdMReM[16] = Foam::Pair<Foam::scalar>(0.495,1164);
CdMReM[17] = Foam::Pair<Foam::scalar>(0.494,1313);
CdMReM[18] = Foam::Pair<Foam::scalar>(0.498,1461);
CdMReM[19] = Foam::Pair<Foam::scalar>(0.503,1613);
CdMReM[20] = Foam::Pair<Foam::scalar>(0.511,1764);
CdMReM[21] = Foam::Pair<Foam::scalar>(0.52,1915);
CdMReM[22] = Foam::Pair<Foam::scalar>(0.529,2066);
CdMReM[23] = Foam::Pair<Foam::scalar>(0.544,2211);
CdMReM[24] = Foam::Pair<Foam::scalar>(0.559,2357);
CdMReM[25] = Foam::Pair<Foam::scalar>(0.575,2500);
CdMReM[26] = Foam::Pair<Foam::scalar>(0.594,2636);
CdMReM[27] = Foam::Pair<Foam::scalar>(0.615,2772);
CdMReM[28] = Foam::Pair<Foam::scalar>(0.635,2905);
CdMReM[29] = Foam::Pair<Foam::scalar>(0.66,3033);
CdMReM[30] = Foam::Pair<Foam::scalar>(0.681,3164);
CdMReM[31] = Foam::Pair<Foam::scalar>(0.7,3293);
CdMReM[32] = Foam::Pair<Foam::scalar>(0.727,3423);
CdMReM[33] = Foam::Pair<Foam::scalar>(0.751,3549);

//Foam::List<Foam::Tuple2<Foam::scalar,Foam::Pair<Foam::scalar> > > allData; // fh, dia, vt
allData.setSize(17);

//allData[0] = Foam::Tuple2<Foam::scalar,Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0003,1.17));
allData[0] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0003,1.17));
allData[1] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0004,1.62));
allData[2] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0005,2.06));
allData[3] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0006,2.47));
allData[4] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0007,2.87));
allData[5] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0008,3.27));
allData[6] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0009,3.67));
allData[7] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.001,4.03));
allData[8] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0012,4.64));
allData[9] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0014,5.17));
allData[10] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0016,5.65));
allData[11] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.0018,6.09));
allData[12] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.002,6.49));
allData[13] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.003,8.06));
allData[14] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.004,8.83));
allData[15] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.005,9.09));
allData[16] = Foam::Tuple2<Foam::scalar, Foam::Pair<Foam::scalar> >(0,Foam::Pair<Foam::scalar>(0.006,9.18));

