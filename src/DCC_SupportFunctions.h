/// Creation Eigen::Sparse_Matrix from file
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols) {

    Eigen::SparseMatrix<double> res(Rows, Cols);
    typedef Eigen::Triplet<double> Tr; // Eigen library class
    std::vector<Tr> tripletList; // Probe vector of triplets

    double i = 0.0, j = 0.0, value = 0.0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
            // cout << i << "\t" << j << "\t" << value << endl;
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong
//Sparse AB matrix
    res.setFromTriplets(tripletList.begin(), tripletList.end());

//Remove all elements and free the memory from the probe vector
    inAN.close();
    tripletList.clear();
    tripletList.shrink_to_fit();

    return res;
}

/// Function for reading triplet list from file
vector<Triplet<double>> TripletsReader(char* SMpath) {
    typedef Eigen::Triplet<double> Tr; // Eigen library class
    std::vector<Tr> tripletList; // Probe vector of triplets

    double i = 0, j = 0, value = 0, t_length = 0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(Tr(i , j, value));
            //   cout << i << "\t" << j << "\t"<< value << "\t" << endl;
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong

    return tripletList;
}

/// Function for reading tuples list from file
vector<vector<double>> VectorVectors4Reader(char* SMpath) {
    vector<vector<double>> res;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            double i = 0.0, j = 0.0, k = 0.0, l = 0;
            inAN >> i >> j >> k >> l;
            res.push_back({i,j,k,l});
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong

    return res;
}

vector<tuple<double, double, double>> TuplesReader(char* SMpath) {
    typedef tuple<double, double, double> Tup; // Eigen library class
    std::vector<Tup> tripletList; // Probe vector of triplets

    double i = 0.0, j = 0.0, value = 0.0, t_length = 0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(make_tuple(i, j, value));
            //   cout << i << "\t" << j << "\t"<< value << "\t" << endl;
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong

    return tripletList;
}
vector<tuple<double, double, double>> dTuplesReader(char* SMpath, unsigned int &size) {
    typedef tuple<double, double, double> Tup; // Eigen library class
    std::vector<Tup> tripletList; // Probe vector of triplets

    double i = 0.0, j = 0.0, value = 0.0, t_length = 0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(make_tuple(i, j, value));
            //   cout << i << "\t" << j << "\t"<< value << "\t" << endl;
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong

    size = tripletList.size();
    return tripletList;
}

/// Creation int std::Vector from file
std::vector<unsigned int> VectorReader(char* FilePath) {
    std::vector<unsigned int> res;
    unsigned int i=0;
    ifstream inCellNumbers(FilePath);
    if (inCellNumbers.is_open()) { //If the file was successfully open, then
        while(inCellNumbers >> i) res.push_back(i);
    } else std:cout << "The file " << FilePath << " cannot be read !" << endl; //If something goes wrong

    return res;
}

/// Creation double std::Vector from file
std::vector<double> dVectorReader(char* FilePath) {
    std::vector<double> res;
    double i = 0.0;
    ifstream inCellNumbers(FilePath);
    if (inCellNumbers.is_open()) { //If the file was successfully open, then
        while(inCellNumbers >> i) res.push_back(i);
    } else std:cout << "The file " << FilePath << " cannot be read" << endl; //If something goes wrong

    return res;
}

/// * The function count the number of Edges possessing with types J0, J1, J2, J3 and J4 for every junction * ///
/// * Calculation Face-Edge index ::                                                                                                     * ///

vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, vector<int> const& TJsTypes) {
    vector<double> res(100,0); /// Up to 100 types of possible TJs types

    for (unsigned int l = 0; l < FES.rows(); l++) // Loop over all Edges
        if (FES.coeff(l,face_number) == 1) res[TJsTypes.at(l)]++;
    // output in the form res[0] = #TJsTypes[0] incident to the face with the number face_number, res[1] = #TJsTypes[1] incident to the face with the number face_number,...

    return res;
}

/// * Function calculates the vector<int> "TJsTypes" of types TJs in the DCC using its FES incidence matrix and special faces sequence (s_faces_sequence) * ///
/// *                                                                                                                                                    * ///
vector<int> EdgesTypesCalc(std::vector<unsigned int> const &CellNumbs, vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &FES)
{
    vector<int> TJsTypes(CellNumbs.at(1),0); // CellNumbs.at(1) is the number of Edges
    for (auto vit: s_faces_sequence) // loop over all Special Faces
        for(int k = 0; k < CellNumbs.at(1); k++) // loop over all Edges
            if (FES.coeff(k, vit) == 1) TJsTypes.at(k)++;

    return TJsTypes;
}

/// DDRX support function :: GFS matrix reading and calculation of new seeds at the centres of GBs
tuple<double, double, double> find_aGBseed(unsigned int facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, vector<tuple<double, double, double>> & AllSeeds_coordinates) {
    tuple <double, double, double> res; // find two grain neighbour for fnumber
    vector<double> xx, yy, zz;
    //Triplet<double> res;     // find two grain neighbour for fnumber

    SpMat GFS(CellNumbs.at(2),CellNumbs.at(3));
    /// GFS matrix reading
    GFS = SMatrixReader(paths.at(6), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains

    vector<int> Grain_neighbours;
    vector<unsigned int> grainIDs;
#pragma omp parallel for // parallel execution by OpenMP
  /*
        grainIDs.clear();
        for (SparseMatrix<double>::InnerIterator it(GFS, facenumb); it; ++it)
            grainIDs.push_back(it.col());
*/
    //normal loop
for (unsigned int j = 0; j < CellNumbs.at(3); ++j) {
        if (GFS.coeff(facenumb, j) == 1) {
            grainIDs.push_back(j);
            if(j<CellNumbs.at(3)-1) {
             for (unsigned int k = j + 1; k < CellNumbs.at(3); ++k)
                if (GFS.coeff(facenumb, k) == 1)
                    grainIDs.push_back(k);
            }
        }
}
//        cout << grainIDs[0] << " " << grainIDs[1] << endl;

    if(grainIDs.size()>0) xx.push_back(get<0>(AllSeeds_coordinates.at(grainIDs[0])));
    if(grainIDs.size()>1) xx.push_back(get<0>(AllSeeds_coordinates.at(grainIDs[1])));
    if(grainIDs.size()>0) yy.push_back(get<1>(AllSeeds_coordinates.at(grainIDs[0])));
    if(grainIDs.size()>1) yy.push_back(get<1>(AllSeeds_coordinates.at(grainIDs[1])));
    if(grainIDs.size()>0) zz.push_back(get<2>(AllSeeds_coordinates.at(grainIDs[0])));
    if(grainIDs.size()>1) zz.push_back(get<2>(AllSeeds_coordinates.at(grainIDs[1])));
    if(grainIDs.size()>1) res = make_tuple(0.5*(xx[0] + xx[1]), 0.5*(yy[0] + yy[1]), 0.5*(zz[0] + zz[1]));
        else if(grainIDs.size()>0) res = make_tuple(xx[0], yy[0], zz[0]);
        else res = make_tuple(0, 0, 0);
//    cout << "facenumb " << facenumb << " " << 0.5*(xx[0] + xx[1]) << "\t" << 0.5*(yy[0] + yy[1]) << "\t"<< 0.5*(zz[0] + zz[1]) << "\t" << endl;
//    cout << "facenumb " << facenumb << " " << get<0>(res) << "\t" << get<1>(res) << "\t"<< get<2>(res) << "\t" << endl;

    return res;

}

bool is_file_exists(const string fileName)
{
    char* charfileName = const_cast<char*>(fileName.c_str());
    std::ifstream infile(charfileName);
    return infile.good();
}

/*
 *     double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0, Face_Entropy_Median = 0 , Face_Entropy_Skrew = 0;

 */

/*
 *     /// Ordinary face sequence
    ordinary_faces_sequence.clear();
    for (auto  itr = State_sVector.begin(); itr != State_sVector.end(); ++itr)
        if (*itr == 0) ordinary_faces_sequence.push_back(distance(State_sVector.begin(), itr));

 */