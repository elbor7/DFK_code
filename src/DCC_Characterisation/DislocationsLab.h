/*kgc=10.0*powf(10,16);
BurgV=2.56*powf(10,-10);
YStress0=1.3*powf(10,8);
ShStress=46*powf(10,9);
Vc=5.8*powf(10,-8);
Kalf=1.6;
RHOD00=powf(10,12);

dRHODm = DeF*(kgc*BurgV*YStress0 + 0.4*ShStress*BurgV*BurgV*kgc*sqrtf(RHOD) - Vc*(RHODm - RHOD00)*sqrtf(RHODim) - Kalf*(2.0*RHODm + RHODim));
dRHODim = DeF*(Vc*(RHODm - RHOD00)*sqrtf(RHODim) - Kalf*RHODim);
//				if (i==0) {				dRHODm = HAGBsFunc[i][0]*(kgc*BurgV*YStress0 + 0.4*ShStress*BurgV*BurgV*kgc*sqrtf(RHOD) - Vc*(RHODm - RHOD00)*sqrtf(RHODim) - Kalf*(2.0*RHODm + RHODim));
//				dRHODim = HAGBsFunc[i][0]*(Vc*(RHODm - RHOD00)*sqrtf(RHODim) - Kalf*RHODim); }
//				else {				dRHODm = (HAGBsFunc[i][0]-HAGBsFunc[i-1][0])*(kgc*BurgV*YStress0 + 0.4*ShStress*BurgV*BurgV*kgc*sqrtf(RHOD) - Vc*(RHODm - RHOD00)*sqrtf(RHODim) - Kalf*(2.0*RHODm + RHODim));
//				dRHODim = (HAGBsFunc[i][0]-HAGBsFunc[i-1][0])*(Vc*(RHODm - RHOD00)*sqrtf(RHODim) - Kalf*RHODim); }

dRHOD = dRHODm + dRHODim;
RHODm += dRHODm;	RHODim += dRHODim;	RHOD += dRHOD;
*/
