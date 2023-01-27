///================================ DCC Face Laplacian module =============================================================///
///=======================================================================================================================///
/** This subfile calculates the combinatorial Laplacian of the special Faces graph with its spectrum **/
///=======================================================================================================================///


vector<double> FaceLaplacian(Eigen::SparseMatrix<double> const& LSF, std::vector<unsigned int> &CellNumbs) {
//REPAIR    cout << LSF.rows() - 2 << endl;
    std::vector<double> Lface_spectrum;
    Eigen::VectorXcd eigval;
    //Eigen::MatrixXcd eigvec;

/// Eigenvalues of the S-Faces laplacian
// Construct matrix operation object using the wrapper class SparseGenMatProd
    SparseGenMatProd<double> op(LSF);
// Construct eigen solver object, requesting all the eigenvalues
    GenEigsSolver<SparseGenMatProd<double>> eigs(op, LSF.rows() - 2, LSF.rows());
    //nev	Number of eigenvalues requested. This should satisfy 1≤𝑛𝑒𝑣≤𝑛−2, where 𝑛 is the size of matrix.
    //ncv	Parameter that controls the convergence speed of the algorithm. Typically a larger ncv means faster convergence, but it may also result in greater memory use and more matrix operations in each iteration. This parameter must satisfy 𝑛𝑒𝑣+2≤𝑛𝑐𝑣≤𝑛, and is advised to take 𝑛𝑐𝑣≥2⋅𝑛𝑒𝑣+1.
// Initialize and compute
    eigs.init();
//    unsigned int nconv = LSF
    eigs.compute(SortRule::LargestMagn);

// Retrieve results
//cout << MatrixXd(LSF) << endl;
    if (eigs.info() == CompInfo::Successful) {
        eigval = eigs.eigenvalues();
//        eigvec = eigs.eigenvectors();
    }

//for (unsigned int cs = 0; cs < (LSF.rows() - 2); ++cs) Lface_spectrum.at(cs) =  eigval.coeffRef(cs,0);
   for (unsigned int i = 0; i < LSF.rows() - 2; ++i) if(real(eigval(i,0)) >= pow(10,-5)) Lface_spectrum.push_back(real(eigval(i,0)));
    for (unsigned int i = 0; i < LSF.rows() - 2; ++i) if(real(eigval(i,0)) < pow(10,-5)) Lface_spectrum.push_back(0.0);

   //cout << Lface_spectrum.size() / (LSF.rows() - 3);
    //Matrix<complex<double>, -1, 1> eigval //    std::cout << "\n Eigenvectors nonZeros :\n" << eigvec.nonZeros() << "\n Eigenvectors size :\n" << eigvec.size()  << std::endl; //                        cout << Laplacian_specialFaces << endl;

    return Lface_spectrum;
} /// End of FaceLaplacian(Laplacian_SFaces)


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * 	//Eigenvalues for DENSE matrices
	MatrixXd LapFEig(Facenumb, Facenumb);
//	MatrixXd LapFEig(LimFaceNumb, LimFaceNumb);

	for (int i = 0; i < Facenumb; i++) for (int j = 0; j < Facenumb; j++) LapFEig(i,j)=LapF[i][j];
	for (int i = 0; i < LimFaceNumb; i++) for (int j = 0; j < LimFaceNumb; j++) LapFEig(i,j)=LapF[i][j];
	EigenSolver<MatrixXd> es(LapFEig);
	cout<<"es(LapFEig) just calculated"<<"    size = "<< es.eigenvalues().size() <<endl;
	MatrixXcd EigF = es.eigenvalues();

	//REsistivity and Conductivity
	RG=0.0;
//	for (int ig = 0; ig < (es.eigenvalues().size()-1); ig++) cout<< real(EigF(ig,0))<<endl;
	for (int ig = 0; ig < (es.eigenvalues().size()-1); ig++) { RG += pow(real(EigF(ig,0)),-1);}
		RG = RG*es.eigenvalues().size();
		SG = pow(RG,-1);

//	ofstream AMFStream;		//	AMFStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)LapF_HAGBs.txt", ios::trunc);
			AMFStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)Conduct.txt", ios::app);
	//		for (int i = 0; i < Facenumb; i++) { for (int j = 0; j < Facenumb; j++) {
	//				if(AMF[i][j]>=0) AMFStream << AMF[i][j] << "\t"; else j++; }			AMFStream << "\n";		}
//			for (int i = 0; i < Facenumb; i++) {for (int j = 0; j < Facenumb; j++) AMFStream << LapF[i][j] << "\t"; AMFStream << "\n"; }
//			for (int i = 0; i < Facenumb; i++) {for (int j = 0; j < Facenumb; j++) AMFStream << LapFEig(i,j) << "\t"; AMFStream << "\n"; }
//			AMFStream << endl << EigF << endl << endl;
			AMFStream << HAGBsFunc[i][1] << "\t" << RG<< "\t" << SG << endl;
			AMFStream.close();	//for(int kii = 0; kii < Facenumb; kii++)	{for(int l = 0; l < 15; l++) if (NewFaceNeighbours[kii][l]>=0) cout<<"   "<<NewFaceNeighbours[kii][l];} cout<<endl;
//	cout<<"EigF file is written"<<endl;  //	system("pause");
 */