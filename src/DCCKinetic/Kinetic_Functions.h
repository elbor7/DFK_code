///================================ A part of the DCC Kinetic module =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new identifications for k-Cells               *
 *  in the DCC Kinetic module. It makes them "special" and takes out of the set of "ordinary" k-Cells.                            **/
///================================================================================================================================///

//int DCC_Kinetic_Wear(double ShearStress, vector<Tup> &Grain_Orientations, Eigen::SparseMatrix<double> const& FES, std::vector<unsigned int> &CellNumbs, char* input_folder, char* output_dir) {
int DCC_Kinetic_Wear(double ShearStress, vector<Tup> &Grain_Orientations, Eigen::SparseMatrix<double> const& FES) {

    /// Input and output directories
    char* odir = const_cast<char*>(output_folder.c_str()); // const_cast for output directory
    char* indir = const_cast<char*>(input_folder.c_str()); // const_cast for output directory

/// Vectors
    /// Normals - reading from file
    string Norm_path = input_folder + "Normal.txt"s;
    char* normals = const_cast<char*>(Norm_path.c_str());
    /// Calculation of the NormalForce for each Face
    std::vector<Tup> NormalForce;
    std::vector<Tup> Norms_tuple = TuplesReader(normals); // Reading from file Normal vectors of all Faces
    for (unsigned long ik = 0; ik < CellNumbs.at(2); ++ik) NormalForce.push_back( { ShearStress*get<0>(Norms_tuple[ik]), ShearStress*get<1>(Norms_tuple[ik]),ShearStress*get<2>(Norms_tuple[ik]) } );
    /// Set grain orientations (Random + From file)
    // I. Random case
    // The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)
    srand((unsigned) time(NULL));

    cout << "Hello there!" << endl;
    return 1;
}

///*============================================================================================*///
///*============================== DCC_Kinetic_Plasticity function =============================*///

//vector<Tup> DCC_Kinetic_Plasticity(Eigen::SparseMatrix<double> const& FES, std::vector<unsigned int> &CellNumbs, char* input_folder, char* output_dir)  {
vector<Tup> DCC_Kinetic_Plasticity( Eigen::SparseMatrix<double> const& FES)  {

    //resultant tuple
    vector<Tup> fraction_stress_temperature;
/// Functions
    vector<unsigned int> Metropolis(vector<vector<double>> &stress_tensor, vector<vector<double>> &norms_vector, vector<vector<double>> &tang_vector, double &Temperature, std::vector<unsigned int> &CellNumbs, long iteration_number, vector<double> &slip_vector, double alpha, double lambda);
    vector<vector<double>> lt_vector(vector<vector<double>> &stress_tensor, vector<vector<double>> norms_vector);

/// Variables
    vector<unsigned int> State_Vector(CellNumbs.at(2),0); // State vector filling with zeros
    /// Other parameters
    vector<vector<double>> tang_vector; // tangential vector to the slip plane (lt)

    /// Normals - reading from file
    vector<vector<double>> norms_vector; // vector of normals
    string Norm_path = input_folder + "normals.txt"s;
    char* normals = const_cast<char*>(Norm_path.c_str());
    // Calculation of the NormalForce for each Face
    std::vector<Tup> Norms_tuple = TuplesReader(normals); // Reading from file Normal vectors of all Faces; TuplesReader() defined in DCC_SupportFunctions.h
    for (unsigned long ik = 0; ik < CellNumbs.at(2); ++ik) norms_vector.push_back( { get<0>(Norms_tuple[ik]), get<1>(Norms_tuple[ik]), get<2>(Norms_tuple[ik]) } );

    /// Slip areas for each Face - reading from file
    vector<double> SAreas, SAreas_m2; // vector of areas : unit + dimension
    string SAreas_path = input_folder + "faces_areas.txt"s;
    char* areas = const_cast<char*>(SAreas_path.c_str());
    SAreas = dVectorReader(areas); // reading from file; dVectorReader() defined in DCC_SupportFunctions.h

    /// Iteration over ShearStresses :: Input parameters
    double ShearStress = 0.0;
    double ShearStress_min = 0.1*pow(10,9), ShearStress_max = 20.0*pow(10,9);
    int simulation_steps = 1000000, dSS = 0;
    /// Iteration over Temperature :: Input parameters
    double Temperature = 293.0;
    /// Material parameters
    double a_latticeCu = 3.628*pow(10,-10); // lattice parameter for bulk COPPER
    double BurgvCu = 2.56*pow(10,-10); // dislocation Burgers vector for COPPER
    double ShearModCu = 46.0*pow(10,9); // Shear modulus for COPPER

    double a_latticeAl = 4.0479*pow(10,-10); // lattice parameter for bulk ALUMINIUM
    double ShearModAl = 26.0*pow(10,9); // Shear modulus for ALUMINIUM
    double BurgvAl = 2.86*pow(10,-10); // Shear modulus for ALUMINIUM

    /// Model parameters
    double Burgv = BurgvCu, ShMod = ShearModCu, a_lattice = a_latticeCu;
    double lambda = 0.0*pow(10,-9)*46.0*pow(10,9), //Shear modulus
             alpha = 1.0*pow(10,-1)*ShMod; ///////
                     //pow(10,-3)*ShMod; /// model coefficient alpha*s^2
    double D_size = 11.0*a_lattice; /// Complex size ///
    vector<double> external_normal = { 0.0, 0.0, 1.0 }; /// normal to the external face
    /// From unit areas to the real ones [m^2]
    for(auto slareas : SAreas) SAreas_m2.push_back(slareas*pow(a_lattice,2));

    ShearStress = ShearStress_min;
    dSS = (ShearStress_max - ShearStress_min)/ (double) simulation_steps;
    cout << "Stress calculation step =\t" << dSS/pow(10,6) << endl;

    for (long i = 0; i < simulation_steps; ++i) {
        ShearStress += dSS;
        /// Stress tensor definition
        double s11 = 0.0, s12 = 0.0, s13 = 0.0, s21 = 0.0, s22 = 0.0, s23 = 0.0, s31 = 0.0, s32 = 0.0, s33 = ShearStress;
        vector<vector<double>> stress_tensor = {{s11, s12, s13}, {s21, s22, s23}, {s31, s32, s33} };
        /// OBTAINING OF NORMAL VECTORS
        tang_vector = lt_vector(stress_tensor, norms_vector);
         //for (auto p : tang_vector) cout << p[0] << "\t" << p[1] << "\t" << p[2] << endl; exit(32);
        /// Slip vectors
        vector<double> slip_vector; // vector of nano-slips mudulus (s)
        for(auto sa : SAreas_m2) slip_vector.push_back(Burgv * sqrt(M_PI/sa)); // nano-slip value s = b * Sqrt(Pi/As)
            //for (auto p : slip_vector) cout << "slip_vector" << SAreas_m2.size() << endl;

        /// METROPOLIS algorithm ///
        /// Iteration number for METROPOLIS algorithm
        long iteration_number = 6.0 * SAreas.size();

        /// =========== Metropolis algorithm ============>>>
        vector<double> smax = *max_element(begin(stress_tensor), end(stress_tensor));
        double s_min = *max_element(begin(smax), end(smax));

        std::fill(State_Vector.begin(), State_Vector.end(), 0);
        if (alpha - s_min < -100.0) std::fill(State_Vector.begin(), State_Vector.end(), 1);
        else if (alpha - s_min > -100.0 && alpha - s_min < 100.0) State_Vector = Metropolis(stress_tensor, norms_vector, tang_vector, Temperature, CellNumbs, iteration_number, slip_vector, alpha, lambda); // Metropolis() function is defined below
       // for(auto state : State_Vector) cout << state << "\t";
       // cout << "End" << endl;
        //cout << "alpha - s_min  =\t" << alpha - s_min << endl;

        /// Plastic strain
        double Plastic_Strain = 0.0; unsigned int state_it = 0;
        for(auto state : State_Vector) {
            //cout << count(begin(State_Vector), end(State_Vector),1) << SAreas_m2.at(state_it) * Burgv * abs(tang_vector.at(state_it)[0]*external_normal[0] + tang_vector.at(state_it)[1]*external_normal[1] + tang_vector.at(state_it)[2]*external_normal[2])/ pow(D_size, 3) << endl;
            if (state == 1) {
                Plastic_Strain += SAreas_m2.at(state_it) * Burgv * abs(tang_vector.at(state_it)[0]*external_normal[0] + tang_vector.at(state_it)[1]*external_normal[1] + tang_vector.at(state_it)[2]*external_normal[2])/ pow(D_size, 3);

            } // end if
            ++state_it;
        } // end of for(state : State_Vector)

        /// Analysis :: slip_fraction and fraction_stress_temperature
        double slip_fraction = std::count(State_Vector.begin(), State_Vector.end(), 1)/ (double) State_Vector.size();
        fraction_stress_temperature.push_back(make_tuple(Plastic_Strain, ShearStress, Temperature));

        //cout << "Stress  =\t" << ShearStress << "\tPlastic strain\t" << Plastic_Strain << endl;
        /// Output and stop!
        if (Plastic_Strain >= 0.002) { cout << "Plastic strain =\t" << Plastic_Strain << "\tYield strength [GPa] =\t" << ShearStress/pow(10,9) << endl; return fraction_stress_temperature;}
      //  cout << "\tSlip fraction =\t" << slip_fraction << "\tPlastic strain =\t" << Plastic_Strain << "\tYield strength [GPa] =\t" << ShearStress/pow(10,9) << endl;

    } // end of for (< calculation_steps)

    cout << "Stress is not in the range!" << endl;
    return fraction_stress_temperature;
}

vector<unsigned int> Metropolis(vector<vector<double>> &stress_tensor, vector<vector<double>> &norms_vector, vector<vector<double>> &tang_vector, double &Temperature, std::vector<unsigned int> &CellNumbs, long iteration_number, vector<double> &slip_vector, double alpha, double lambda){
    /// I. Constant initial state initialisation
    //zero vectors required for Processing_Random() input
    std::vector <unsigned int> SlipState_Vector(CellNumbs.at(2),0), s_faces_sequence(CellNumbs.at(2),0);
    double Rc = 8.31; //gas constant

    // ================> Initial p = 0.5 Face seeds
    Processing_Random( SlipState_Vector, s_faces_sequence, 0.5);

    // 1. Loop over all the Faces
    srand((unsigned) time(NULL)); // The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)

    for (long i = 0; i < iteration_number; ++i) {
    // 2. Random choice of a new Face
    long NewSlipNumber = rand() % (CellNumbs.at(2) -2); // Random generation of the boundary number in the range from 0 to CellNumbs.at(2)

    // 3. If dH < 0 energetically favourable -> Accept trial and change the type
        if (SlipState_Vector.at(NewSlipNumber) == 1) {
            SlipState_Vector.at(NewSlipNumber) = 0;
        } //if
     // 4. Else if dH > 0 consider the acceptance probability
        else if (SlipState_Vector.at(NewSlipNumber) == 0) {
    // Model variables
            double Snt = 0.0, s2 = 0.0;
            vector<vector<double>> sik{ {0,0,0}, {0,0,0}, {0,0,0} };
//            cout << NewSlipNumber << "\t" << slip_vector.at(NewSlipNumber) << endl;
            //exit(457);
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    sik[i][j] = 0.5 * slip_vector.at(NewSlipNumber) * (norms_vector.at(NewSlipNumber)[i] * tang_vector.at(NewSlipNumber)[j] +
                                                norms_vector.at(NewSlipNumber)[j] * tang_vector.at(NewSlipNumber)[i]);
                    Snt += stress_tensor[i][j] * sik[i][j];

                    s2 += 0.5 * pow(slip_vector.at(NewSlipNumber),2) * (norms_vector.at(NewSlipNumber)[i] * tang_vector.at(NewSlipNumber)[j]*norms_vector.at(NewSlipNumber)[i] * tang_vector.at(NewSlipNumber)[j] +
                            norms_vector.at(NewSlipNumber)[i] * tang_vector.at(NewSlipNumber)[j]*norms_vector.at(NewSlipNumber)[j] * tang_vector.at(NewSlipNumber)[i]);
                }
            }
            /// New ACCEPTANCE PROBABILITY
            double P_ac = 0.0;
            if (s2 != 0) {
                double P_ac = exp(-(alpha * pow(slip_vector.at(NewSlipNumber),2) - Snt - lambda * slip_vector.at(NewSlipNumber)) / (Rc * Temperature));
                if (P_ac > 1) P_ac = 1.0;
            } else P_ac = 0;
           // cout << "Alpha\t" << alpha * pow(slip_vector.at(NewSlipNumber),2) - Snt  << "\t SNT\t" << Snt << "\tProb\t" << exp(-(alpha * s2 - Snt - lambda * slip_vector.at(NewSlipNumber)) / (Rc * Temperature)) << endl;

            double rv = (rand() / (RAND_MAX + 1.0)); // Generate random value in the range [0,1]
                if (rv <= P_ac) {
    //   if (P_ac > 0) cout << "P_ac =\t" << P_ac << "\trandom number\t" << rv << endl;
                    SlipState_Vector.at(NewSlipNumber) = 1;
                } else SlipState_Vector.at(NewSlipNumber) = 0;

        } // end of  else if (SlipState_Vector.at(NewSlipNumber) == 0)

    } // for loop (i < iteration_number)

    return SlipState_Vector;
} // end of Metropolis

vector<vector<double>> lt_vector(vector<vector<double>> &stress_tensor, vector<vector<double>> norms_vector) {
    vector<vector<double>> tv, tnv, ttv, lt; // traction vector tv and its normal (tnv) and tangent (ttv) components + tangential vector to the slip plane lt

//for (auto kl : stress_tensor.at(1))    cout << kl << endl;
    for (unsigned int fnumb = 0; fnumb < norms_vector.size(); ++fnumb) { // loop over all the slip elements in the complex

        double tv0 = inner_product(stress_tensor.at(0).begin(), stress_tensor.at(0).end(),
                                   norms_vector.at(fnumb).begin(), 0);
        double tv1 = inner_product(stress_tensor.at(1).begin(), stress_tensor.at(1).end(),
                                   norms_vector.at(fnumb).begin(), 0);
        double tv2 = inner_product(stress_tensor.at(2).begin(), stress_tensor.at(2).end(),
                                   norms_vector.at(fnumb).begin(), 0);
        tv.push_back({tv0, tv1, tv2});

        double tnv0 = inner_product(tv.at(fnumb).begin(), tv.at(fnumb).end(), norms_vector.at(fnumb).begin(), 0) *
                      norms_vector[fnumb][0];
        double tnv1 = inner_product(tv.at(fnumb).begin(), tv.at(fnumb).end(), norms_vector.at(fnumb).begin(), 0) *
                      norms_vector[fnumb][1];
        double tnv2 = inner_product(tv.at(fnumb).begin(), tv.at(fnumb).end(), norms_vector.at(fnumb).begin(), 0) *
                      norms_vector[fnumb][2];
        tnv.push_back({tnv0, tnv1, tnv2});
        double ttv0 = tv0 - tnv0, ttv1 = tv1 - tnv1, ttv2 = tv2 - tnv2;
        ttv.push_back({ttv0, ttv1, ttv2});

        double  lt0 = 0.0, lt1 = 0.0, lt2 = 0.0;
        if (inner_product(ttv.at(fnumb).begin(), ttv.at(fnumb).end(), ttv.at(fnumb).begin(), 0.0L) != 0) {
            lt0 =
                    ttv0 / sqrt(inner_product(ttv.at(fnumb).begin(), ttv.at(fnumb).end(), ttv.at(fnumb).begin(), 0.0L));
             lt1 =
                    ttv1 / sqrt(inner_product(ttv.at(fnumb).begin(), ttv.at(fnumb).end(), ttv.at(fnumb).begin(), 0.0L));
             lt2 =
                    ttv2 / sqrt(inner_product(ttv.at(fnumb).begin(), ttv.at(fnumb).end(), ttv.at(fnumb).begin(), 0.0L));
        }
        lt.push_back({lt0, lt1, lt2});

    } // end of for(fnumb < norms_vector.size())

    return lt;
} /// end of lt_vector() function

//std::vector <unsigned int> DCC_Kinetic_cracking(std::vector<unsigned int> &s_faces_sequence, std::vector<unsigned int> &CellNumbs, Eigen::SparseMatrix<double> const& AFS, Eigen::SparseMatrix<double> const& FES, std::vector<char*> paths, char* input_folder, char* output_dir) {
std::vector <unsigned int> DCC_Kinetic_cracking(std::vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const& AFS, Eigen::SparseMatrix<double> const& FES) {

    std::vector <unsigned int> crack_faces_sequence, S_crackVector(CellNumbs.at(2), 0); // sequence of the cracked Faces
    vector<int> TJsTypes(CellNumbs.at(1), 0), TJsCrackTypes(CellNumbs.at(1), 0);
    vector<double> newCrack_neigh_TJs, newCrack_neigh_Faces; // only for neighbouring faces
    std::vector <double> Face_weight(CellNumbs.at(2),0), Face_rGO_index(CellNumbs.at(2),0), Face_crack_index(CellNumbs.at(2),0); // indexes for all GBs
    std::vector <double> lsc_rGO_energy(CellNumbs.at(2),0), lsc_crack_energy(CellNumbs.at(2),0); // energies related with the local stress concentrators
    std::vector <double> Face_energy(CellNumbs.at(2),0), Face_current_energy(CellNumbs.at(2),0); // energies for all GBs
    double crack_fraction = 0.0;
    double sface_energy_rGO = 0.5*pow(10,-3), sface_energy_matrix = 1.0*pow(10,-3);
    double  lsc_rGO = 4.0*sface_energy_rGO, lsc_crack = 2.0*sface_energy_rGO; // energy values related with the stress concentrators
    double kB = 1.3807*pow(10,-23); // Boltzmann constant

    /// Initial energies
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i) Face_energy.at(i) = sface_energy_matrix; // assignments of the matrix surface energies
    for (auto ks : s_faces_sequence) Face_energy.at(ks) = sface_energy_rGO; // assignments of the rGO defect surface energies

    //Calculation of the TJs types
    TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
    /// GB_indices calculation
    for (unsigned int fn = 0; fn < CellNumbs.at(2); ++fn) {
        vector<double> j_types_neigh_fractions = GBIndex(fn, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
        //REPAIR   for(auto kl : j_types_neigh_fractions) cout << " " <<kl ; cout << endl;

        // rGO index calculation
        Face_rGO_index.at(fn) = (j_types_neigh_fractions.at(0) + 2.0 * j_types_neigh_fractions.at(1) + 3.0 * j_types_neigh_fractions.at(2)) / 3.0;

    } // end of for ( fn < CellNumbs.at(2))

    /// Local stress concentrators and related energies
    // Face weights calculation
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i)
        for (int l = 0; l < AFS.rows(); l++) // Loop over all Faces
            if (AFS.coeff( l, i) == 1) Face_weight.at(i)++;

    // Energies
    for (unsigned int j = 0; j < CellNumbs.at(2); ++j) lsc_rGO_energy.at(j) = lsc_rGO*(Face_rGO_index.at(j)/Face_weight.at(j));
    for (unsigned int k = 0; k < CellNumbs.at(2); ++k) lsc_crack_energy.at(k) = lsc_crack*(Face_crack_index.at(k)/Face_weight.at(k));

    /// Final surface energy
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i) Face_energy.at(i) = Face_energy.at(i) -= lsc_rGO_energy.at(i);

    /// Crack face sequence creation
    /// (1) FAST case : cracks do NOT affect the process of the further fracture`
    //     for (unsigned int i = 0; i < CellNumbs.at(2); ++i){ ;}

    /// (2) SLOW case : cracks affect the process of the further fracture`
    double NewCrackNumb = 0;
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i) { // crack loop over all the GBs
        /// New surface energy
        // Local stress concentrators and related energies
        for (unsigned int k = 0; k < CellNumbs.at(2); ++k) lsc_crack_energy.at(k) = lsc_crack*(Face_crack_index.at(k)/Face_weight.at(k));
        // Current surface energy
        for (unsigned int l = 0; l < CellNumbs.at(2); ++l) Face_current_energy.at(l) = Face_energy.at(l) - lsc_crack_energy.at(l);

        /// Number of element
        NewCrackNumb = std::min_element(std::begin(Face_energy), std::end(Face_energy)) - std::begin(Face_energy); // gives index of the max element
        if (find(crack_faces_sequence.begin(), crack_faces_sequence.end(), NewCrackNumb) == crack_faces_sequence.end()) {
            crack_faces_sequence.push_back(NewCrackNumb); // if the element is not in the sequence - add
            S_crackVector.at(NewCrackNumb) = 1;
            Face_energy.at(NewCrackNumb) = 100000.0; /// arbitrarily large value !!!
        } // end if

        /// Neighbours of the converted Face
        newCrack_neigh_TJs.clear();
        newCrack_neigh_Faces.clear();
        for (int k = 0; k < CellNumbs.at(1); ++k) // Loop over all the Edges
            if (FES.coeff(k, NewCrackNumb) == 1) newCrack_neigh_TJs.push_back(k);
        if (newCrack_neigh_TJs.size() > 0) for (auto itr: newCrack_neigh_TJs) TJsCrackTypes.at(itr)++;

        for (int m = 0; m < CellNumbs.at(2); ++m) // Loop over all the Faces
            if (AFS.coeff(NewCrackNumb, m) == 1 && S_crackVector.at(m) == 0) newCrack_neigh_Faces.push_back(m);

        if (newCrack_neigh_Faces.size() > 0) {
            for (auto fn: newCrack_neigh_Faces) {
                vector<double> j_types_crack_n_fractions = GBIndex(fn, FES,TJsCrackTypes); //Types (up to 100 kinds) of the edges incident to the considered Face

                /// New crack index
                Face_crack_index.at(fn) = (j_types_crack_n_fractions.at(0) + 2.0 * j_types_crack_n_fractions.at(1) +
                                           3.0 * j_types_crack_n_fractions.at(2)) / 3.0;
            } // end for (auto fn: new_crack_neigh_Faces)
        } // end if (new_crack_neigh_Faces.size() > 0)

        crack_fraction = crack_faces_sequence.size() / (double) CellNumbs.at(2);
//REPAIR cout << crack_fraction << endl;
    } // end for ( i < CellNumbs.at(2) )

    return crack_faces_sequence;
}