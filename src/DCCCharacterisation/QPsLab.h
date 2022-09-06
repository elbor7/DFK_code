/*//-------------------------------------------------------------------------------------
//Ji-Jj analyser
//--------------------------------------------------------------------------------------
				for (int in = 0; in < Edgenumb; in++)
				for (int j = 0; j < 100; j++)
					JEdgeNeigh[in][j] = -10;

for(int ih = 0; ih < Edgenumb; ih++)	{
		knn=0;
//Сначала мы выясняем тип самой грани
	J1count=0; 	JN0Count=0;
	for (int lk = 0; lk < Facenumb; lk++) { if(MFE1[lk][ih] == 2) JN0Count++; if(MFE1[lk][ih] >= 1) J1count++; }
		if(J1count > 1) JEdgeNeigh[ih][0] = JN0Count;
//			{cout<<"  il= "<<ih<<"  JEdgeNeigh0= "<<JEdgeNeigh[ih][0]<<endl; system("pause");}
		//Сначала в 0-ячейку собственный тип ребра            //
				}

for(int ih = 0; ih < Edgenumb; ih++)	{
		knn=1;
//Затем проходим всех ее соседей
	for (int lmn = 0; lmn < 100; lmn++) {
		if(EdgeNeighbours[ih][lmn]>=0)  knumb = EdgeNeighbours[ih][lmn];
			else knumb=-1;
//------------------------------------------------------
		J1count=0; 		JNCount=0;
		if(knumb>=0) { for (int lk = 0; lk < Facenumb; lk++) { if(MFE1[lk][knumb] == 2) JNCount++; if(MFE1[lk][knumb] >= 1) J1count++; }
//-----------------------------------------------------
		if(J1count > 1) JEdgeNeigh[ih][knn++] = JNCount;
			else JEdgeNeigh[ih][knn+1]=-1; //knn++

		}			}		}

//Вывод в файл JEdgeNeigh.txt
if(i==10) {
	JEdgeN.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)JEdgeNeigh10.txt", ios::trunc);
			for (int il = 0; il < Edgenumb; il++) {
				for (int jl = 0; jl < 100; jl++) {
					if(JEdgeNeigh[il][jl]>=0) JEdgeN << JEdgeNeigh[il][jl] << "\t";
				}
			JEdgeN <<"Edge number="<<il<<"\n";
		}
			JEdgeN.close();
			}

	for (int ip = 0; ip < 10; ip++)
		for (int j = 0; j < 10; j++)
				JEN[ip][j] = 0;

//Анализ матрицы JEdgeNeigh[i][j]
for (int il = 0; il < Edgenumb; il++) //для каждой грани
	for (int nl = 0; nl < 10; nl++) //мы перебираем все варианты какой она может быть
		 if(JEdgeNeigh[il][0] == nl)  for (int mj = 1; mj < 100; mj++) //НАЧАЛО С 1, ЧТОБЫ НЕ УЧИТЫВАТЬ САМУ ГРАНЬ КАК СОСЕДА// и если она оказалась определенного типа, то мы перебираем всех ее соседей
										for (int nk = 0; nk < 10; nk++) if(JEdgeNeigh[il][mj] == nk) JEN[nl][nk]++; //так что если сосед оказывается также определенного типа, то мы заносим их связь в матрицу JEN


//Полсчет "мощности" каждого тройного стыка с учетом соседей
for (int il = 0; il < Edgenumb; il++) {
	TJpow[il]=0; NEneigh=0; for (int jl = 0; jl < 100; jl++)	if(JEdgeNeigh[il][jl]>=0) {TJpow[il]+= JEdgeNeigh[il][jl]; NEneigh++; };
	TJpow[il] = TJpow[il]/NEneigh;
//Только ненулевые стыки
	TJpow2[il]=0; NEneigh2=0; for (int jl = 0; jl < 100; jl++)	if(JEdgeNeigh[il][jl]>=0 && JEdgeNeigh[il][0]>0) {TJpow2[il]+= JEdgeNeigh[il][jl]; NEneigh2++; };
	TJpow2[il] = TJpow2[il]/NEneigh2;
}


//Подсчет среднего по комплексу и дисперсии DISPERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NEneighSum=0.0; NEneighAv=0; for (int il = 0; il < Edgenumb; il++) if(TJpow[il]>=0) {NEneighSum += TJpow[il];  NEneighAv++;}
NEneighSum = NEneighSum/NEneighAv;

//Только ненулевые стыки
NEneighSum2=0.0; NEneighAv2=0; for (int il = 0; il < Edgenumb; il++) if(TJpow2[il]>0) {NEneighSum2 += TJpow2[il];  NEneighAv2++;}
NEneighSum2 = NEneighSum2/NEneighAv2;

//Dispersion
SNEneighSum=0.0; NEneighAv=0; for (int il = 0; il < Edgenumb; il++) if(TJpow[il]>=0) {SNEneighSum += (TJpow[il]*TJpow[il]); NEneighAv++;}
SNEneighSum = SNEneighSum/(NEneighAv);
//Только ненулевые стыки
SNEneighSum2=0.0; NEneighAv2=0; for (int il = 0; il < Edgenumb; il++) if(TJpow2[il]>0) {SNEneighSum2 += (TJpow2[il]*TJpow2[il]); NEneighAv2++;}
SNEneighSum2 = SNEneighSum2/(NEneighAv2);

if((SNEneighSum - NEneighSum*NEneighSum) >0) DispNES1 = sqrtl(SNEneighSum - NEneighSum*NEneighSum);
else DispNES1 = sqrtl(-SNEneighSum +NEneighSum*NEneighSum);
if((SNEneighSum2 - NEneighSum2*NEneighSum2) >0) DispNES2 = sqrtl(SNEneighSum2 - NEneighSum2*NEneighSum2);
else  DispNES2 = sqrtl(-SNEneighSum2 +NEneighSum2*NEneighSum2);

	JEdgeN.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)TJspow.txt", ios::app);
		JEdgeN << HAGBsFunc[i][0]<<"\t"<<25.0*powf(10,6)/(sqrtf(RHOD))<<"\t"<<30.0*powf(10,6)/(HAGBsFunc[i][1]*sqrtf(RHOD))<<"\t"<< 30.0*powf(10,6)*NEneighSum/(NEneighSum2*HAGBsFunc[i][1]*sqrtf(RHOD))<<"\t"<< 3.0*27.0*powf(10,6)/(NEneighSum2*HAGBsFunc[i][1]*sqrtf(RHOD));

//				if(NEneighSum2>0)  JEdgeN <<HAGBsFunc[i][0]<<"\t"<<NEneighSum<<"\t"<<NEneighSum2<<"\t"<<DispNES1<<"\t"<<DispNES2<<"\t"<<2.2*powl(NEneighSum,(-1.0/1.0))<<"\t"<<2.2*powl(NEneighSum2,(-1.0/1.0)); //<<"\t"<<SNEneighSum<<"\t"<<SNEneighSum2;
//				else JEdgeN <<HAGBsFunc[i][0]<<"\t"<<NEneighSum<<"\t"<<0.0;
//				 JEdgeN <<"Accumulated Strain =  "<<HAGBsFunc[i][0]<<"\n";
//				 JEdgeN<< NEneighSum << "\t";
//Мощность всех узлов
//				for(int il = 0; il < Edgenumb; il++) if(TJpow[il]>=0) JEdgeN<<TJpow[il] << "\t";
			JEdgeN <<"\n";
			JEdgeN.close();

//--------------------------------------------------------------------------------------

//Вывод в файл JEN.txt
		JENStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)JEN.txt", ios::app);
			JENStream <<AcStrain<<endl;
			for (int i = 0; i < 10; i++) {
				for (int j = 0; j < 10; j++) if(JEN[i][j]>=0)  JENStream << JEN[i][j] << "\t";
						//			JENStream << JEN[i][j] << "\t";
									JENStream <<"\n";
			}
// Очевидно, что при таком алгоритме диагональные элементы учитываются дважды, а также считаем число всех элементов и делим потом на него, то есть
			SumJEN=0; Jenii=0; Jenij=0;
			for (int ic = 0; ic < 10; ic++)
				for (int jc = 0; jc < 10; jc++) {SumJEN += 0.5*JEN[ic][jc]; if(ic==jc) Jenii+=0.5*JEN[ic][jc]; else Jenij+=0.5*JEN[ic][jc]; }
//cout<< SumJEN <<endl;

			JENStream <<"\n"<<"\n"<<"\n";
		JENStream.close();
//***cout<<"(HAGBs)JEN.txt has been created"<<endl;

		JENStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)Jlk.txt", ios::app);
//				if(i==0)	JENStream << "STRAIN" << "\t"<<"ii"<<"\t"<<"ij"<<"\t"<<"00"<<"\t"<<"11"<<"\t"<<"22"<<"\t"<<"33"<<"\t"<<"01"<<"\t"<<"02"<<"\t"<<"03"<<"\t"<<"12"<<"\t"<<"13"<<"\t"<<"23"<<endl;
				if(i==0)	JENStream << "STRAIN" << "\t"<<"ii"<<"\t"<<"ij"<<"\t"<<endl;
				JENStream << HAGBsFunc[i][0]<<"\t"<<Jenii*100.0/SumJEN<<"\t"<< Jenij*100.0/SumJEN;
//				JENStream << HAGBsFunc[i][0]<<"\t"<<JEN[0][0] + JEN[1][1] + JEN[2][2] + JEN[3][3] + JEN[4][4]<<"\t"<< JEN[0][1] + JEN[0][2] + JEN[0][3] + JEN[1][2] + JEN[1][3] + JEN[2][3] + JEN[2][4] + JEN[3][4] << "\t"<< JEN[0][0] << "\t"<< JEN[1][1] << "\t"<< JEN[2][2] << "\t"<< JEN[3][3] << "\t"<< JEN[0][1] << "\t"<< JEN[0][2] << "\t"<< JEN[0][3] << "\t"<< JEN[1][2] << "\t"<< JEN[1][3] << "\t"<< JEN[2][3];     //			 JENStream << HAGBsFunc[i][0]<<"\t"<<JEN[0][0] + JEN[1][1] + JEN[2][2] + JEN[3][3]<<"\t"<< JEN[0][1] + JEN[1][0] + JEN[0][2] + JEN[2][0] + JEN[0][3] + JEN[3][0] + JEN[1][2] + JEN[2][1] + JEN[1][3] + JEN[3][1] + JEN[2][3] + JEN[3][2]<< "\t"<< JEN[0][1] + JEN[1][0]<< "\t"<< JEN[0][2] + JEN[2][0]<< "\t"<< JEN[0][3] + JEN[0][3]<< "\t"<< JEN[1][2] + JEN[2][1]<< "\t"<< JEN[1][3] + JEN[3][1]<< "\t"<< JEN[2][3] + JEN[3][2];
//			 JENStream << HAGBsFunc[i][0]<<"\t"<<JEN[0][0] + JEN[1][1] + JEN[2][2] + JEN[3][3]<<"\t"<< JEN[0][1] + JEN[0][2] + JEN[0][3] + JEN[1][2] + JEN[1][3] + JEN[2][3] << "\t"<< JEN[0][1] << "\t"<< JEN[0][2] << "\t"<< JEN[0][3] << "\t"<< JEN[1][2] << "\t"<< JEN[1][3] << "\t"<< JEN[2][3]; 	// 			 JENStream << 0*JEN[0][0] + 2*JEN[1][1] + 4*JEN[2][2] + 6*JEN[3][3]<<"\t"<< JEN[0][1] + JEN[1][0] + 2*JEN[0][2] + 2*JEN[2][0] + 3*JEN[0][3] + 3*JEN[3][0] + 3*JEN[1][2] + 3*JEN[2][1] + 4*JEN[1][3] + 4*JEN[3][1] + 5*JEN[2][3] + 5*JEN[3][2]<< "\t"<< JEN[0][1] + JEN[1][0]<< "\t"<< 2*JEN[0][2] + 2*JEN[2][0]<< "\t"<< 3*JEN[0][3] + 3*JEN[0][3]<< "\t"<< 3*JEN[1][2] + 3*JEN[2][1]<< "\t"<< 4*JEN[1][3] + 4*JEN[3][1]<< "\t"<< 5*JEN[2][3] + 5*JEN[3][2];
					JENStream <<"\n";
		JENStream.close();
*/
