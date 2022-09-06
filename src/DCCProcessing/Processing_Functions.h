///================================ A part of the DCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new identifications for k-Cells               *
 *  in the DCC Processing module. It makes them "special" and takes out of the set of "ordinary" k-Cells.                            **/
///================================================================================================================================///

/// (1) Totally random choice of element with # New2CellNumb from the list of numbers {0 to SCellsNumb}
unsigned int NewCellNumb_R(unsigned int OCellsNumb){ // Random generation of a 2-Cell number
    unsigned int New2CellNumb;
    New2CellNumb = rand() % (OCellsNumb-1); // Random generation of the boundary number in the range from 0 to OrdinaryCellNumbs.size()-1
    return New2CellNumb;
}

/// (2) The Random generation process function
/*!
 * @param S_Vector
 * @param s_faces_sequence
 * @param max_sFaces_fraction
 */
int Processing_Random( std::vector<unsigned int> &S_Vector,  std::vector<unsigned int> &s_faces_sequence, double max_sFaces_fraction) {
///================================================================= 'R' =======================================================================////
/// ====================================================>  Random generation process  <========================================================////
///===========================================================================================================================================////

///=============================================================================================================================================////
/// =====> Initial initialisation with the previous calculation step (if any) based on the _special_faces_sequence_

// OrdinaryCellNumbs is just a tricky way to acceleration the random process
    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(2), 1); // Vector of the size equal to the total number of faces in DCC initialised with '0's
    for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++) OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces
    for( unsigned int itr : S_Vector) if(itr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + itr); // !!! Delete its element from the vector decreasing its size BUT

    // unsigned int OCellAmount = std::count(SpecialCellNumbs->begin(), SpecialCellNumbs->end(), 0);
    double ordinary_faces_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(2);
    double special_faces_fraction = 1.0 - ordinary_faces_fraction;
    if (special_faces_fraction >= max_sFaces_fraction) return 0;

/// ================= Loop before max_sFaces_fraction fractions of special cells  =======================>
    /// The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)

    do { // do{ ... }while(output_step) loop starting point
        srand((unsigned) time(NULL));

        int New2CellNumb = 0;
        New2CellNumb = NewCellNumb_R(OrdinaryCellNumbs.size());
//REPAIR cout << OrdinaryCellNumbs.size() << "\t\t" << New2CellNumb << "\t\t" << OrdinaryCellNumbs.at(New2CellNumb) << endl;
        /// Random generation of types !!! with IDs < number_of_types
        int NewFaceType = 1;
        if (number_of_types > 1) NewFaceType = rand() % number_of_types; // Random chose for the chosen boundary to be assigned over all special types
        /// Changes in vectors from Main function
        S_Vector.at(OrdinaryCellNumbs.at(New2CellNumb)) = NewFaceType;
        s_faces_sequence.push_back(OrdinaryCellNumbs.at(New2CellNumb));

        // It is a key tricky point for the fast Random generation process: vector decreases after each turn from an ordinary to special face BUT we chose all the boundary one by one because their initial numeration stored in the SpecialCellNumbs[] vector !
        OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + New2CellNumb); // !!! Delete its element from the vector decreasing its size

        // Special and Ordinary Faces fraction calculation
        ordinary_faces_fraction = OrdinaryCellNumbs.size() / (double) CellNumbs.at(2);
        special_faces_fraction = 1.0 - ordinary_faces_fraction;

    } while(special_faces_fraction < max_sFaces_fraction); /// End of the Random generation process
/// Closing and deleting
    // Remove all elements anf free the memory from the probe OrdinaryCellNumbs vector
    OrdinaryCellNumbs.clear();
    OrdinaryCellNumbs.shrink_to_fit();

    return 0;
} // end  of Random

/// (3) Maximum entropy generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
//int Processing_maxEntropy(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, double max_sFaces_fraction, int number_of_types, std::vector<unsigned int> &CellNumbs, std::vector<char*> paths) {
int Processing_maxEntropy(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence) {
///=============================================================================================================================================////
///========================================================================= 'S' =================================================================////
/// ====================================================>  Maximum entropy production process   <========================================================////
///=============================================================================================================================================////

    /// Local functions declarations
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Sstr0 = 0, Sstr00 = 0, SstrN = 0, Sstr0N = 0, Configuration_Face_Entropy = 0, Face_Entropy_Median = 0, Face_Entropy_Skrew = 0;

    // Obtaining Faces (coloumns) - Edges (rows) incidence matrix from file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1), CellNumbs.at(2)); // Edges-Faces
    SpMat AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose())); // Full matrix instead of triagonal

    double special_faces_fraction = s_faces_sequence.size() / (double) CellNumbs.at(2);
    double ordinary_faces_fraction = 1.0 - special_faces_fraction;
    if (special_faces_fraction >= max_sFaces_fraction) {
        cout << "(!) Early exit: Fraction of special faces is already >= MAX values from config.txt" << endl;
        return 0;
    }
    else if (special_faces_fraction < 0.05) {
        /// ================> Initial Face seeds (initial state for the MAX entropy production algorithm)
        Processing_Random(S_Vector, s_faces_sequence, 0.05);
///***function
    }

    /// ========== Calculation configuration entropy at each calculation step ===============>
    /// Vectors for Edges types and Edges-related configuration entropy
    vector<int> TJsTypes(CellNumbs.at(1), 0); // vector int in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values
    // Zeroing coefficients of TJsTypes vector
    std::fill(TJsTypes.begin(), TJsTypes.end(), 0);
    /// TJs type calculations
    TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
///***function
//REPAIR for(auto kl : TJsTypes) cout << " " <<kl ; cout << endl; exit(10);
    J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    Configuration_Face_Entropy = 0;

    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
    J0 = CellNumbs.at(1) - J1 - J2 - J3;
    Jall = CellNumbs.at(1);
// Conversion from numbers to fractions | (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    j0 = J0 / Jall;
    j1 = J1 / Jall;
    j2 = J2 / Jall;
    j3 = J3 / Jall;
    /// Configuration Entropy related with Faces
    Configuration_Face_Entropy = -(j0 * log2(j0) + j1 * log2(j1) + j2 * log2(j2) + j3 * log2(j3));

    /// ======== Loop over all Faces ===========>
    vector<double> EntropyIncreaseList(CellNumbs.at(2), 0); // vector with values of configuration entropy increases at conversion of each Face
    for (unsigned int k = 0; k < CellNumbs.at(2); k++) { // loop over all Faces in DCC
        double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_Entropy = 0, CFace_EntropyIncrease = 0;

        if (S_Vector.at(k) == 0) { // Loop over each still ORDINARY element neighbours
            J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;

            vector<double> j_types_neigh_fractions = GBIndex(k, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
///***function
//REPAIR   for(auto kl : j_types_neigh_fractions) cout << " " <<kl ; cout << endl;
            // Values before conversion
            J00 = j_types_neigh_fractions.at(0);
            J10 = j_types_neigh_fractions.at(1);
            J20 = j_types_neigh_fractions.at(2);
//REPAIR    cout << " J00= " << J00<< " J10= " << J10 << " J20= " << J20 << endl;
            // Values after conversion
            J1N = J00;
            J2N = J10;
            J3N = J20;
            // The entropy increase calculation for a given Face
            CFace_EntropyIncrease = (J0 * log2(J0 + pow(10, -30)) - J00 * log2(J00 + pow(10, -30)))
                                    + (J1 * log2(J1 + pow(10, -30)) - J10 * log2(J10 + pow(10, -30)) +
                                       J1N * log2(J1N + pow(10, -30)))
                                    + (J2 * log2(J2 + pow(10, -30)) - J20 * log2(J20 + pow(10, -30)) +
                                       J2N * log2(J2N + pow(10, -30)))
                                    + (J3 * log2(J3 + pow(10, -30)) + J3N * log2(J3N + pow(10,-30)));
//REPAIR  cout  << "\t\t" <<  CFace_EntropyIncrease << "\t\t" << endl;

            // The result of one iteration
            EntropyIncreaseList.at(k) = CFace_EntropyIncrease;

        } // if OrdinaryCells (S_Vector.at(Face) == 0)
    } // for (..k < CellNumbs.at(2)..)

    double New2CellNumb = 0; // Only one possible Face type (binary model)
    /// Number of element giving the maximum increase in configuration entropy at its conversion
  New2CellNumb = std::max_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
///min  New2CellNumb = std::min_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
//REPAIR        cout << s_faces_sequence.size() << "   " << New2CellNumb << endl;

    // Then all the corresponding maps chain
    S_Vector.at((unsigned int) New2CellNumb) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

    /// Add the new element to s_faces_sequence if it is still not here
    if (find(s_faces_sequence.begin(), s_faces_sequence.end(), New2CellNumb) == s_faces_sequence.end())
        s_faces_sequence.push_back(New2CellNumb);

    EntropyIncreaseList.at(New2CellNumb) = 0.0; ///working_max

///=============================================================================================================================================////
/// ================= S_MAX LOOP over the fractions of special cells [0.05,1] =======================>
    vector<double> new_neigh_TJs, new_neigh_Faces;
    do { // do{ ... }while(output_step) loop starting point

        /// Neighbours of the converted Face
        new_neigh_TJs.clear();
        new_neigh_Faces.clear();
        for (int k = 0; k < CellNumbs.at(1); ++k) // Loop over all the edges
            if (FES.coeff(k, New2CellNumb) == 1) new_neigh_TJs.push_back(k);
        if (new_neigh_TJs.size() > 0) for (auto itr: new_neigh_TJs) TJsTypes.at(itr)++;

        for (int m = 0; m < CellNumbs.at(2); ++m) // Loop over all the Faces
            if (AFS.coeff(New2CellNumb, m) == 1 && S_Vector.at(m) == 0) new_neigh_Faces.push_back(m);

        if (new_neigh_Faces.size() > 0) {
            for (auto itr: new_neigh_Faces) {
                vector<double> j_types_neigh_fractions = GBIndex(itr, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face

                /// New EntropyIncreaseList changes
                double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;
                // Values before conversion
                J00 = j_types_neigh_fractions.at(0);
                J10 = j_types_neigh_fractions.at(1);
                J20 = j_types_neigh_fractions.at(2);
                // Values after conversion
                J1N = J00;
                J2N = J10;
                J3N = J20;
                // The entropy increase calculation for a given Face
                CFace_EntropyIncrease = (J0 * log2(J0 + pow(10, -30)) - J00 * log2(J00 + pow(10, -30)))
                                        + (J1 * log2(J1 + pow(10, -30)) - J10 * log2(J10 + pow(10, -30)) + J1N * log2(J1N + pow(10, -30)))
                                        + (J2 * log2(J2 + pow(10, -30)) - J20 * log2(J20 + pow(10, -30)) + J2N * log2(J2N + pow(10, -30)))
                                        + (J3 * log2(J3 + pow(10, -30)) + J3N * log2(J3N + pow(10,-30)));
                /// The result of the one iteration
                EntropyIncreaseList.at(itr) = CFace_EntropyIncrease; // for(auto kl :EntropyIncreaseList) cout << kl << endl;

            } // for (auto itr : new_neigh_Faces )
        } //if

        double New2CellNumb = 0; // Only one possible Face type (binary model)
        /// Number of element giving the maximum increase in configuration entropy at its conversion
        New2CellNumb = std::max_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element

///min         unsigned int imp = 0; for( auto str : EntropyIncreaseList) { if( str == 0) EntropyIncreaseList.at(imp) = pow(10,10); ++imp; }
///min        New2CellNumb = std::min_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
//REPAIR     cout  << "\t\t" <<  New2CellNumb << "\t\t" << EntropyIncreaseList.at(New2CellNumb) << endl;

        // Then all the corresponding maps chain
        S_Vector.at((unsigned int) New2CellNumb) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

        /// Add the new element to s_faces_sequence if it is still not here
        if(find(s_faces_sequence.begin(),s_faces_sequence.end(),New2CellNumb) == s_faces_sequence.end())
            s_faces_sequence.push_back(New2CellNumb);

        EntropyIncreaseList.at(New2CellNumb) = 0.0; ///working_max
///min        EntropyIncreaseList.at(New2CellNumb) = pow(10,10);
//REPAIRS       cout << New2CellNumb << "  " << EntropyIncreaseList.at(New2CellNumb) << endl; //for(auto kl :s_faces_sequence) cout << " " <<kl ;
//INFO cout << "\tFace fraction = \t" << s_faces_sequence.size() / (double) CellNumbs.at(2) << endl; //exit(10);

        /// Special and Ordinary Faces fraction calculation
        special_faces_fraction = s_faces_sequence.size() / (double) CellNumbs.at(2);
        ordinary_faces_fraction = 1.0 - special_faces_fraction;

    } while(special_faces_fraction < max_sFaces_fraction); /// End of the Random generation process
//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;

/// Closing and deleting
     TJsTypes.clear();
     TJsTypes.shrink_to_fit();
    EntropyIncreaseList.clear();
    EntropyIncreaseList.shrink_to_fit();
    new_neigh_TJs.clear();
    new_neigh_TJs.shrink_to_fit();
    new_neigh_Faces.clear();
    new_neigh_Faces.shrink_to_fit();
    FES.makeCompressed();
    AFS.makeCompressed();

    return 0;
} // end of S_max
//REPAIRS/////////        for(auto kl :S_Vector) cout << " " <<kl ; cout << endl; //exit(10);

/// (4) i(p) index govern processed
int Processing_ipIndex(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, int index_type, double ip_index) {
///=============================================================================================================================================////
///========================================================================= 'X' (ip INDEX) =================================================================////
/// ====================================================>  i(p) Index-based production process   <========================================================////
///=============================================================================================================================================////
    /// Local functions declarations
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double l20 = 0, l21 = 0, l22 = 0, l23 = 0, l200 = 0, l210 = 0, l21N = 0, l220 = 0, l22N = 0, l230 = 0, l23N = 0;
    double Sstr0 = 0, Sstr00 = 0, SstrN = 0, Sstr0N = 0, Configuration_Face_Entropy = 0, Face_Entropy_Median = 0, Face_Entropy_Skrew = 0;

    // Obtaining Faces (coloumns) - Edges (rows) incidence matrix from file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1), CellNumbs.at(2)); // Edges-Faces
    SpMat AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose())); // Full matrix instead of triagonal

    double special_faces_fraction = s_faces_sequence.size() / (double) CellNumbs.at(2);
    double ordinary_faces_fraction = 1.0 - special_faces_fraction;
    if (special_faces_fraction >= max_sFaces_fraction) {
        cout << "(!) Early exit: Fraction of special faces is already >= MAX values from config.txt" << endl;
        return 0;
    }
    else if (special_faces_fraction < 0.05) {
        /// ================> Initial Face seeds (initial state for the MAX entropy production algorithm)
        Processing_Random(S_Vector, s_faces_sequence, 0.05);
///***function
    }

    /// ========== Calculation configuration entropy at each calculation step ===============>
    /// Vectors for Edges types and Edges-related configuration entropy
    vector<int> TJsTypes(CellNumbs.at(1), 0); // vector int in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values
    // Zeroing coefficients of TJsTypes vector
    std::fill(TJsTypes.begin(), TJsTypes.end(), 0);
    /// TJs type calculations
    TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
///***function
//REPAIR for(auto kl : TJsTypes) cout << " " <<kl ; cout << endl; exit(10);
    J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    Configuration_Face_Entropy = 0;

    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
    J0 = CellNumbs.at(1) - J1 - J2 - J3;
    Jall = CellNumbs.at(1);
// Conversion from numbers to fractions | (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    j0 = J0 / Jall;
    j1 = J1 / Jall;
    j2 = J2 / Jall;
    j3 = J3 / Jall;
    /// Configuration Entropy related with Faces
    Configuration_Face_Entropy = -(j0 * log2(j0) + j1 * log2(j1) + j2 * log2(j2) + j3 * log2(j3));

    /// ======== Loop over all Faces ===========>
    vector<double> EntropyIncreaseList(CellNumbs.at(2), 0); // vector with values of configuration entropy increases at conversion of each Face
    for (unsigned int k = 0; k < CellNumbs.at(2); k++) { // loop over all Faces in DCC
        double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_Entropy = 0, CFace_EntropyIncrease = 0;

        if (S_Vector.at(k) == 0) { // Loop over each still ORDINARY element neighbours
            J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;

            vector<double> j_types_neigh_fractions = GBIndex(k, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
///***function
//REPAIR   for(auto kl : j_types_neigh_fractions) cout << " " <<kl ; cout << endl;
            // Values before conversion
            J00 = j_types_neigh_fractions.at(0)/ Jall;
            J10 = j_types_neigh_fractions.at(1)/ Jall;
            J20 = j_types_neigh_fractions.at(2)/ Jall;
//REPAIR    cout << " J00= " << J00<< " J10= " << J10 << " J20= " << J20 << endl;
            // Values after conversion
            J1N = J00;
            J2N = J10;
            J3N = J20;
            if (j0 != 0) l20 = log2(j0); else l20 = 0;   if (J00 != 0) l200 = log2(J00); else l200 = 0;
            if (j1 != 0) l21 = log2(j1); else l21 = 0;   if (J10 != 0) l210 = log2(J10); else l210 = 0; if (J1N != 0) l21N = log2(J1N); else l21N = 0;
            if (j2 != 0) l22 = log2(j2); else l22 = 0;   if (J20 != 0) l220 = log2(J20); else l220 = 0; if (J2N != 0) l22N = log2(J2N); else l22N = 0;
            if (j3 != 0) l23 = log2(j3); else l23 = 0;   if (J30 != 0) l230 = log2(J30); else l230 = 0; if (J3N != 0) l23N = log2(J3N); else l23N = 0;

            // The entropy increase calculation for a given Face
            CFace_EntropyIncrease = (j0 * l20 - J00 * l200) + (j1 * l21 - J10 * l210 + J1N * l21N) + (j2 * l22 - J20 * l220 + J2N * l22N) + (j3 * l23 + J3N *l23N);
            /*            CFace_EntropyIncrease = (j0 * log2(j0 + pow(10, -30)) - J00 * log2(J00 + pow(10, -30)))
                                    + (j1 * log2(j1 + pow(10, -30)) - J10 * log2(J10 + pow(10, -30)) +
                                       J1N * log2(J1N + pow(10, -30)))
                                    + (j2 * log2(j2 + pow(10, -30)) - J20 * log2(J20 + pow(10, -30)) +
                                       J2N * log2(J2N + pow(10, -30)))
                                    + (j3 * log2(j3 + pow(10, -30)) + J3N * log2(J3N + pow(10,-30)));
            CFace_EntropyIncrease = (J0 * log2(J0 + pow(10, -30)) - J00 * log2(J00 + pow(10, -30)))
                                    + (J1 * log2(J1 + pow(10, -30)) - J10 * log2(J10 + pow(10, -30)) +
                                       J1N * log2(J1N + pow(10, -30)))
                                    + (J2 * log2(J2 + pow(10, -30)) - J20 * log2(J20 + pow(10, -30)) +
                                       J2N * log2(J2N + pow(10, -30)))
                                    + (J3 * log2(J3 + pow(10, -30)) + J3N * log2(J3N + pow(10,-30)));
//REPAIR  cout  << "\t\t" <<  CFace_EntropyIncrease << "\t\t" << endl;
*/
            // The result of one iteration
            EntropyIncreaseList.at(k) = - CFace_EntropyIncrease;

        } // if OrdinaryCells (S_Vector.at(Face) == 0)
    } // for (..k < CellNumbs.at(2)..)

    double New2CellNumb = 0; // Only one possible Face type (binary model)
    /// Number of element giving the maximum increase in configuration entropy at its conversion
    New2CellNumb = std::min_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
///min  New2CellNumb = std::min_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
//REPAIR        cout << s_faces_sequence.size() << "   " << New2CellNumb << endl;

    // Then all the corresponding maps chain
    S_Vector.at((unsigned int) New2CellNumb) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

    /// Add the new element to s_faces_sequence if it is still not here
    if (find(s_faces_sequence.begin(), s_faces_sequence.end(), New2CellNumb) == s_faces_sequence.end())
        s_faces_sequence.push_back(New2CellNumb);

    EntropyIncreaseList.at((unsigned int)New2CellNumb) = 0;

///=============================================================================================================================================////
/// ================= S_LOOP over the fractions of special cells [0.05,1] =======================>
    vector<double> new_neigh_TJs, new_neigh_Faces;
    double CNb = New2CellNumb;
    do { // do{ ... }while(output_step) loop starting point

        /// Neighbours of the converted Face
        new_neigh_TJs.clear();
        new_neigh_Faces.clear();
        for (int k = 0; k < CellNumbs.at(1); ++k) // Loop over all the edges
            if (FES.coeff(k, CNb) == 1) new_neigh_TJs.push_back(k);

        if (new_neigh_TJs.size() > 0) for(auto itr: new_neigh_TJs) TJsTypes.at(itr)++;

        for (int m = 0; m < CellNumbs.at(2); ++m) // Loop over all the Faces
            if (AFS.coeff(CNb, m) == 1 && S_Vector.at(m) == 0) new_neigh_Faces.push_back(m);

        if (new_neigh_Faces.size() > 0) {
            for (auto itr: new_neigh_Faces) {
                vector<double> j_types_neigh_fractions;
                j_types_neigh_fractions = GBIndex(itr, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face

                J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
                Configuration_Face_Entropy = 0;

                J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
                J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
                J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
                J0 = CellNumbs.at(1) - J1 - J2 - J3;
                Jall = (double) CellNumbs.at(1);

// Conversion from numbers to fractions | (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
                j0 = J0 / Jall;
                j1 = J1 / Jall;
                j2 = J2 / Jall;
                j3 = J3 / Jall;

                /// New EntropyIncreaseList changes
                double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;
                // Values before conversion
                J00 = j_types_neigh_fractions.at(0)/ Jall;
                J10 = j_types_neigh_fractions.at(1)/ Jall;
                J20 = j_types_neigh_fractions.at(2)/ Jall;
                // Values after conversion
                J1N = J00;
                J2N = J10;
                J3N = J20;
                // The entropy increase calculation for a given Face
                if (j0 != 0) l20 = log2(j0); else l20 = 0;   if (J00 != 0) l200 = log2(J00); else l200 = 0;
                if (j1 != 0) l21 = log2(j1); else l21 = 0;   if (J10 != 0) l210 = log2(J10); else l210 = 0; if (J1N != 0) l21N = log2(J1N); else l21N = 0;
                if (j2 != 0) l22 = log2(j2); else l22 = 0;   if (J20 != 0) l220 = log2(J20); else l220 = 0; if (J2N != 0) l22N = log2(J2N); else l22N = 0;
                if (j3 != 0) l23 = log2(j3); else l23 = 0;   if (J30 != 0) l230 = log2(J30); else l230 = 0; if (J3N != 0) l23N = log2(J3N); else l23N = 0;
                CFace_EntropyIncrease = (j0 * l20 - J00 * l200) + (j1 * l21 - J10 * l210 + J1N * l21N) + (j2 * l22 - J20 * l220 + J2N * l22N) + (j3 * l23 + J3N *l23N);

/*
                CFace_EntropyIncrease = (j0 * log2(j0 + pow(10, -30)) - J00 * log2(J00 + pow(10, -30)))
                                        + (j1 * log2(j1 + pow(10, -30)) - J10 * log2(J10 + pow(10, -30)) + J1N * log2(J1N + pow(10, -30)))
                                        + (j2 * log2(j2 + pow(10, -30)) - J20 * log2(J20 + pow(10, -30)) + J2N * log2(J2N + pow(10, -30)))
                                        + (j3 * log2(j3 + pow(10, -30)) + J3N * log2(J3N + pow(10,-30)));
 */
                /// The result of the one iteration
                EntropyIncreaseList.at(itr) = - CFace_EntropyIncrease; // for(auto kl :EntropyIncreaseList) cout << kl << endl;

            } // for (auto itr : new_neigh_Faces )
        } //if

        double New2CellNumb = 0; // Only one possible Face type (binary model)
        double min_val = 10;
        /// Number of element giving the maximum increase in configuration entropy at its conversion
unsigned int ipl = 0;
// REPAIR for ( auto pl : new_neigh_TJs ) cout << pl << "  "; cout << "new list : "<< endl;
for ( auto pl : EntropyIncreaseList ) {
    if (pl < min_val && pl > 0) {
        min_val = pl;
        New2CellNumb = (double) ipl;
    }
    ++ipl;
}

//REPAIR     cout  << "\t\t" <<  New2CellNumb << "\t\t" << EntropyIncreaseList.at(New2CellNumb) << endl;

        // Then all the corresponding maps chain
        S_Vector.at((unsigned int) New2CellNumb) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

        /// Add the new element to s_faces_sequence if it is still not here
        if(find(s_faces_sequence.begin(),s_faces_sequence.end(),New2CellNumb) == s_faces_sequence.end())
            s_faces_sequence.push_back(New2CellNumb);

 // REPAIR       cout << New2CellNumb << "  " << EntropyIncreaseList.at((unsigned int)New2CellNumb) << endl; //for(auto kl :s_faces_sequence) cout << " " <<kl ;
        EntropyIncreaseList.at((unsigned int) New2CellNumb) = 0;
///min        EntropyIncreaseList.at(New2CellNumb) = pow(10,10);
//REPAIRS
//INFO cout << "\tFace fraction = \t" << s_faces_sequence.size() / (double) CellNumbs.at(2) << endl; //exit(10);
        /// Special and Ordinary Faces fraction calculation
        special_faces_fraction = s_faces_sequence.size() / (double) CellNumbs.at(2);
        ordinary_faces_fraction = 1.0 - special_faces_fraction;
        CNb = New2CellNumb; //New2CellNumb!!!

    } while(special_faces_fraction < max_sFaces_fraction); /// End of the Random generation process
//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;

/// Closing and deleting
    TJsTypes.clear();
    TJsTypes.shrink_to_fit();
    EntropyIncreaseList.clear();
    EntropyIncreaseList.shrink_to_fit();
    new_neigh_TJs.clear();
    new_neigh_TJs.shrink_to_fit();
    new_neigh_Faces.clear();
    new_neigh_Faces.shrink_to_fit();
    FES.makeCompressed();
    AFS.makeCompressed();

    return 0;
} /// end of Processing_ipIndex()

/// (5) DDRX process
/*!
 *
 * @param State_Vector
 * @param special_faces_sequence
 * @return
 */
//int Processing_DDRX(std::vector<unsigned int>  &State_Vector, std::vector<unsigned int>  &special_faces_sequence, double max_sFaces_fraction, std::vector<char*> const paths, int number_of_types, std::vector<unsigned int> &CellNumbs) {
int Processing_DDRX(std::vector<unsigned int>  &State_Vector, std::vector<unsigned int>  &special_faces_sequence) {
///=============================================================================================================================================////
///========================================================================= 'D' =================================================================////
/// ==================================================================>  DDRX process   <========================================================////
///=============================================================================================================================================////
    /// Function
    tuple<double, double, double> find_aGBseed(unsigned int Facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, vector<tuple<double, double, double>> & AllSeeds_coordinates);
    /// Streams
    ofstream NewSeedsStream;
    NewSeedsStream.open(paths.at(8), ios::trunc);
    NewSeedsStream << "New generated seeds\t" << endl;
    NewSeedsStream.close();

    /// Vector
//    vector<Triplet<double>> AllSeeds_coordinates, NewSeeds_coordinates;
    vector<tuple<double, double, double>> AllSeeds_coordinates;
    tuple<double, double, double> NewSeed_coordinates = make_tuple(0, 0, 0);

    /// Face Seeds reading from file Seeds.txt
    AllSeeds_coordinates = TuplesReader(paths.at(7));

/// ============== Constants and model parameters ===================
    double kb = 1.23* pow(10,-23);
    double GSav = 200.0*pow(10,-9), Temp = 300, Stress = 4.0*pow(10.0,6), Rcr = 0.5*GSav/10.0,
            Q_energy = 213.0*pow(10,1), R_const = 8.31;
/// Conversion probability
    double prob_seed = exp(-Q_energy/(R_const*Temp));
    //    double prob_seed = exp(- Stress*(13.0*pow(Rcr,3)/3.0)/(kb*Temp));
    cout << "DDRX module:: prob_seed =\t" << prob_seed << endl;

    srand((unsigned) time(NULL)); // seed for random
    // Output stream opening
    NewSeedsStream.open(paths.at(8), ios::app);

    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber) {
        double rv = (rand() / (RAND_MAX + 1.0));
        //   cout << rv << "\t" << prob_seed << endl;
        if (rv <= prob_seed) {
            NewSeed_coordinates = find_aGBseed(fnumber, paths, CellNumbs, AllSeeds_coordinates);
            NewSeedsStream << fnumber << "\t" << get<0>(NewSeed_coordinates) << "\t" << get<1>(NewSeed_coordinates) << "\t" << get<2>(NewSeed_coordinates) << endl;
        }
    }
    NewSeedsStream.close();

    return 0;
}

std::vector <unsigned int> Smax_sequence_reader(char* SFS_dir) {
    std::vector <unsigned int> s_faces_sequence;
    cout << "-----------------------------------------------------------------------------------------------------------------------------"s << endl;
    cout << "Warning!!: special_faces_sequence successfully loaded from file:\t"s << SFS_dir << " because the Processing is OFF"s << endl;
    cout << "------------------------------------------------------------------------------------------------------------------------------"s << endl;
    s_faces_sequence = VectorReader(SFS_dir); //all Faces

    return s_faces_sequence;
// REPAIR   cout << "S_size =\t" << special_faces_sequence.size() << endl;

}

/// Heap
/*
 * int Processing_ExperimentalData(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, double max_sFaces_fraction, int number_of_types, std::vector<unsigned int> &CellNumbs, std::vector<char*> paths) {
    if ( ProcessingON(confpath, time_step_one) && *Processing_type == 'E') {
        string  pass1path = input_folder + "pass_1_misorientation.txt"s; char* p1path = const_cast<char*>(pass1path.c_str());
        string  pass1_B2path = input_folder + "pass_1_B2.txt"s; char* p1B2path = const_cast<char*>(pass1_B2path.c_str());
        string  pass2path = input_folder + "pass_2_misorientation.txt"s; char* p2path = const_cast<char*>(pass2path.c_str());
        string  pass2_B2path = input_folder + "pass_2_B2.txt"s; char* p2B2path = const_cast<char*>(pass2_B2path.c_str());
        string  pass4path = input_folder + "pass_4_misorientation.txt"s; char* p4path = const_cast<char*>(pass4path.c_str());
        string  pass4_B2path = input_folder + "pass_4_B2.txt"s; char* p4B2path = const_cast<char*>(pass4_B2path.c_str());
        string  pass8path = input_folder + "pass_8_misorientation.txt"s; char* p8path = const_cast<char*>(pass8path.c_str());
        string  pass8_B2path = input_folder + "pass_8_B2.txt"s; char* p8B2path = const_cast<char*>(pass8_B2path.c_str());

        vector<*char> e_paths = {p1path, p1B2path, p2path, p2B2path, p4path, p4B2path, p8path, p8B2path};

    }
    return 0;
}
*/