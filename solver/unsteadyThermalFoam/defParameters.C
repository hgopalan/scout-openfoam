Foam::List<Foam::Pair<scalar> > CdMReM; // first is Cd and second is Re
CdMReM.setSize(34);

CdMReM[0] = Foam::Pair<scalar>(15,1.8);
CdMReM[1] = Pair<scalar>(4.2,9.61);
CdMReM[2] = Pair<scalar>(2.4,23.4);
CdMReM[3] = Pair<scalar>(1.66,43.2);
CdMReM[4] = Pair<scalar>(1.28,68.7);
CdMReM[5] = Pair<scalar>(1.07,98.9);
CdMReM[6] = Pair<scalar>(0.926,134);
CdMReM[7] = Pair<scalar>(0.815,175);
CdMReM[8] = Pair<scalar>(0.729,220);
CdMReM[9] = Pair<scalar>(0.671,269);
CdMReM[10] = Pair<scalar>(0.607,372);
CdMReM[11] = Pair<scalar>(0.57,483);
CdMReM[12] = Pair<scalar>(0.545,603);
CdMReM[13] = Pair<scalar>(0.528,731);
CdMReM[14] = Pair<scalar>(0.517,866);
CdMReM[15] = Pair<scalar>(0.504,1013);
CdMReM[16] = Pair<scalar>(0.495,1164);
CdMReM[17] = Pair<scalar>(0.494,1313);
CdMReM[18] = Pair<scalar>(0.498,1461);
CdMReM[19] = Pair<scalar>(0.503,1613);
CdMReM[20] = Pair<scalar>(0.511,1764);
CdMReM[21] = Pair<scalar>(0.52,1915);
CdMReM[22] = Pair<scalar>(0.529,2066);
CdMReM[23] = Pair<scalar>(0.544,2211);
CdMReM[24] = Pair<scalar>(0.559,2357);
CdMReM[25] = Pair<scalar>(0.575,2500);
CdMReM[26] = Pair<scalar>(0.594,2636);
CdMReM[27] = Pair<scalar>(0.615,2772);
CdMReM[28] = Pair<scalar>(0.635,2905);
CdMReM[29] = Pair<scalar>(0.66,3033);
CdMReM[30] = Pair<scalar>(0.681,3164);
CdMReM[31] = Pair<scalar>(0.7,3293);
CdMReM[32] = Pair<scalar>(0.727,3423);
CdMReM[33] = Pair<scalar>(0.751,3549);

List<Tuple2<scalar,Pair<scalar> > allData; // fh, dia, vt
allData.setSize(17);

allData[0] = Tuple2(0,Pair<scalar>(0.0003,1.17));
