///******************************************************************************************************************************///
///*                                                                                                                            *///
///******************************************************************************************************************************///
///********************************** Dr Elijah Borodin, Manchester, UK, Spring 2022   *****************************************///
///*****************************************************************************************************************************///
///****************************************   DCC Design tool (c)   ************************************************************///
///*****************************************************************************************************************************///

/// Standard C++ libraries (STL):
///------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <tuple>
#include <numeric>
#include <algorithm>
#include <thread>
//#include <execution>

///-------------------------------
/// Attached user defined C++ libraries (must be copied in the directory for STL):
///------------------------------------------
// #include <Eigen/Core>
// #include <Eigen/Dense>
#include <Eigen/SparseCore> // Eigen source: https://eigen.tuxfamily.org/ (2022)
#include <Spectra/MatOp/SparseGenMatProd.h> // Spectra source: https://spectralib.org/ (2022)
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
///------------------------------------------
using namespace std; //Standard namespace
using namespace Eigen; //Eigen namespace
using namespace Spectra; //Spectra namespace

/// Triplets in the form T = T(i, j, value), where i and j are element's a(i, j) indices in the corresponding dense matrix and the third variable is its value
typedef Triplet<double> Tr; // <Eigen library class> Declares a triplet's type with name - Tr
typedef SparseMatrix<double> SpMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles with name - SpMat

/// Declaration of GLOBAL variables (can be used in all the project modules and libraries)
int dim; // problem dimension (2D or 3D) and the number of types of "special" elemets in the considering problem
double design_number = 0;
std::vector<unsigned int> CellNumbs; // the vector named CellNumbs with the number of k-cells of different types; will be read from file
vector<char*> paths; // the vector with pathes to input and output directories and file
double max_sFaces_fraction, max_cFaces_fraction; // maximum fractions of special and fractured (cracked) faces (2-cells)
string input_folder, output_folder; // input and output directories readed from the file config.txt file
std::vector<vector<unsigned int>> special_face_design; // is a list of special_faces_sequence as an output of the module. The first line here is always the random case (zero-design), then following Smax (1-sequence) and Smin (2-sequence), and then all the designs in between Smax and Smin
std::vector<unsigned int> face_strip_distribution; // a vector containing length distribution of special faces (strips)

/// PCC Complex Geometry ::
/// Coordinates
vector<tuple<double, double, double>> vertex_coordinates_vector, face_coordinates_vector, grain_coordinates_vector;

//std::vector<tuple<double, double, double>> vertex_coordinates, grain_barycenter_coordinates, face_barycenter_coordinates;
// Complex size: (lengths) Lx, Ly, Lz
double Lx_size = 1.0, Ly_size = 1.0, Lz_size = 1.0; // initial values, then they can be read from file
//double Lx_size = pow(10.0,-6.0), Ly_size = pow(10.0,-6.0), Lz_size = pow(10.0,-6.0); // initial values, then they can be read from file

/// Measures
std::vector<double> grain_volumes_vector, face_areas_vector, edge_lengths_vector;
std::vector<double> sample_size_vector = {1.0, 1.0, 1.0}; // L_x [m], L_y [m],L_z [m]

/// PCC Complex Energies ::
std::vector <double> face_elastic_energies;
std::vector<double> GB_SE_vector, GB_EEE_vector, GB_CIE_vector, GB_BLE_vector, GB_CLE_vector;

/// All global OFFSTREAMS
ofstream OutFLfile, OutJFile, OutJ2File, OutSFile, OutS2File, OutPowersADistributions, OutAvLengthsADistributions, OutAgglStatistics, OutElCondfile, OutCrackEnergies_file;
//(-)vector<bool> SChar_config; // Characterisation module configuration

/// Declaration of FUNCTIONS, please see the function bodies at the end of the main file
std::vector<double> confCount(char* config, string &Subcomplex_type, string &Processing_type, string &Kinetic_type, string &input_folder, string &output_folder); // Read and output the initial configuration from the config.txt file
string SIMULATION_MODE(char* config); // Check the SIMULATION MODE type in the config.txt file
bool SubcomplexON(char* config, bool time_step_one); // Check the Subcomplex (Section) module status (On/Off) in the config.txt file
bool ProcessingON(char* config, bool time_step_one); // Check the Processing module status (On/Off) in the config.txt file
bool KineticON(char* config, bool time_step_one); // Check the Kinetic module status (On/Off) in the config.txt file
bool CharacterisationON(char* config, bool time_step_one); // Check the Structure Characterisation module status (On/Off) in the config.txt file
bool WriterON(char* config, bool time_step_one); // Check the Structure Writer module status (On/Off) in the config.txt file
void eraseSubStr(std::string & mainStr, const std::string & toErase); // support (technical) function: erase the first occurrence of given substring from main string

/// Including all the main project libraries
/* Various useful functions (must be here - first in this list! ) */
#include "DCC_SupportFunctions.h"

/* Objects library contain classes of various objects related to PCC substructures and elements */
#include "DCC_Objects/DCC_Objects.h"

/* Processing module assigned special id's for the various DCC elements (Faces, Edges,...) */
/* Output: module generates the list of sequences of special element's id's corresponding to different generation principles (random, maximum entropy, etc.). */
/* Based on each sequences in the list, the State vector of special elements (2-cells) can be generated */
#include "DCC_Processing/DCC_Processing.h"

/* Kinetic module assign some new values for the scalar or vector variables defined on the DCC elements (Faces, Edges,...) or new types identificators (IDs) of k-cells */
/* As its output module generates one or several sequences of "fractured" element's id's; the new "fractured" State vector of Faces (2-cells) can be generated */
#include "DCC_Kinetic/DCC_Kinetic.h"

/* Multiphysics module assign various physical quantities (energies, temperature, electrical conductivity) to the k-Cells of PCC subcomplexes */
#include "DCC_Multiphysics/DCC_Multiphysics.h"

/* Sections module calculates reduced PCC subcomplexes (including plain cuts) inheriting (reduced) special element sequences and state vectors of the original (processed) PCC */
#include "DCC_Section/DCC_Subcomplex.h"

/* Characterisation module calculates various structural characteristics of special substructures defined on the DCC elements */
#include "DCC_Characterisation/StructureCharacterisation.h"

/* Writer module perform formatted output of various data structures generated by other modules */
#include "DCC_Writer/DCC_Writer.h"

/// PCC Special Elements
// #1# Macro-cracks
std::vector <macrocrack> large_cracks_vector;

///* ........................................................................................    Main    ................................................................ *///
int main() {
    cout << "-------------------------------------------------------------------------------------" << endl;

/// File path to the configuration profile
// Still not working in Windows:(
    using std::__fs::filesystem::current_path; // to obtain current working directory
    string MainPath = current_path(); // MainPath variables
    eraseSubStr(MainPath, "cmake-build-debug"s); // to delete .../cmake-build-debug/ from the path
    string config = MainPath + "config.txt"s; char* confpath = const_cast<char*>(config.c_str()); // file path with config.txt
    bool time_step_one = 1; // (technical support variable) id for the first step of iteration MAX fraction of Faces

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
    // The source directory and simulation type from file config.txt
    string S_type; // 'P' or 'H' :: This char define the subsection type: 'P' for the whole Plane cut, 'H' for the half-plane cut like a crack
    string P_type; // 'R' or 'S' :: This char define the process type: 'R' for Random, 'S' for maximum configuration Entropy production /// It the plans 'D', 'I' or 'E' :: 'D' for DDRX process (DCC retessellations with the new seeds), 'I' for the Index-based generation modes, 'E' for experimental data obtained by EBSD
    string K_type; // 'F' :: This char define the Kinetic process type: 'F' for the Ising-like model of Fracture /// In the plans 'W' or 'P' :: 'W' for the 3D one-layer film, 'P' for the Ising-like model of Plasticity

/// ConfigVector contains all the control variables of the program readed from the config.txt
    vector<double> ConfigVector = confCount(confpath, S_type, P_type, K_type, input_folder, output_folder);
    char* indir = const_cast<char*>(input_folder.c_str()); // const_cast for input directory //    char* odir = const_cast<char*>(output_folder.c_str()); // const_cast for output directory
    dim = ConfigVector.at(0); // space dimension of the problem (dim = 2 or 3);
    design_number = ConfigVector.at(1);
    std::vector<unsigned int> kface_sequence; // special Faces sequence for the Kinetic module
/// Below the file names with the sparse DCC matrices which must already exit in the input_folder and have the same names (!)
/// They can be obtain by the DCC Generator tool (https://github.com/PRISBteam/Voronoi_DCC_Analyser) based on the Neper output (*.tess file) (https://neper.info/) or DCC Structure Generator tool (https://github.com/PRISBteam/DCC_Structure_Generator) for plasticity problems
    string ssd0 = input_folder + "A0.txt"s, ssd1 = input_folder + "A1.txt"s, ssd2 = input_folder + "A2.txt"s, ssd3 = input_folder + "A3.txt"s, ssd4 = input_folder + "B1.txt"s, ssd5 = input_folder + "B2.txt"s, ssd6 = input_folder + "B3.txt"s,
            seeds = input_folder + "seeds.txt"s; //, NewSeeds = input_folder + "NewSeeds/NewSeeds.txt"s;
//The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
    paths.push_back(const_cast<char*>(ssd0.c_str())); paths.push_back(const_cast<char*>(ssd1.c_str())); paths.push_back(const_cast<char*>(ssd2.c_str())); paths.push_back(const_cast<char*>(ssd3.c_str()));
    paths.push_back(const_cast<char*>(ssd4.c_str())); paths.push_back(const_cast<char*>(ssd5.c_str())); paths.push_back(const_cast<char*>(ssd6.c_str()));
    paths.push_back(const_cast<char*>(seeds.c_str()));  //paths.push_back(const_cast<char*>(NewSeeds.c_str()));

    // Checking file names
    char* number_of_cells; string  ncells;
    // File path with the amounts of the different cells  (1st line for Nodes, 2nd line for Edges, 3rd line for Faces and (in 3D) 4th line for grains)
    if (is_file_exists(input_folder + "voro_Ncells.txt"s)) { ncells = input_folder + "voro_Ncells.txt"s; number_of_cells = const_cast<char*>(ncells.c_str()); }
    else if (is_file_exists(input_folder + "delau_Ncells.txt"s)) { ncells = input_folder + "delau_Ncells.txt"s; number_of_cells = const_cast<char*>(ncells.c_str()); }
    else { ncells = input_folder + "number_of_cells.txt"s; number_of_cells = const_cast<char*>(ncells.c_str()); }

    /// Useful ofstreams cleaning for data output
    OutJFile.open(output_folder + "TJp_random_theory.txt"s, ios::trunc); OutJFile.close();
    OutJ2File.open(output_folder + "TJp_crystalline_theory.txt"s, ios::trunc); OutJ2File.close();
    OutSFile.open(output_folder + "Sp_random_theory.txt"s, ios::trunc); OutSFile.close();
    OutJ2File.open(output_folder + "Sp_crystalline_theory.txt"s, ios::trunc); OutS2File.close();
    OutAgglStatistics.open(output_folder + "Stat_Agglomerations.txt"s, ios::trunc); OutAgglStatistics.close(); // Agglomerations statistics output

/// Output paths.vector to console out
    int npath = 0;
    cout << "_____________________________________________________________________________________" << endl;
    for (auto m : paths) cout <<"[" << npath++ << "]" << " paths:\t" << m << endl;

/// Principal variables
///_____________________________________________________________________________________

/// Read all the number of cells from file number_of_cells
// :: CellNumbs vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Grains number
    CellNumbs = VectorReader(number_of_cells); // VectorReader is a function from the DCC_SupportFunctions.h library; "number_of_cells" here if the PATH to file

    /// All grain coordinates
    string GCpath_string = input_folder + "grain_seeds.txt"s;
    char* GCpath = const_cast<char*>(GCpath_string.c_str());
    vector<tuple<double, double, double>> grain_coordinates(CellNumbs.at(3));
    grain_coordinates = TuplesReader(GCpath);
    /// Set grain_vector_coordinates
    for(auto grain_coord_tuple : grain_coordinates)
        grain_coordinates_vector.push_back(grain_coord_tuple);

//REPAIR    cout << "Size of grain_seeds: " << gs_size << " CellNumbs.at(3): " << CellNumbs.at(3) << endl;

    /// All face coordinates /// LONG CALCULATION HERE
    string FSpath_string = input_folder + "face_seeds.txt"s;
    char* FSpath = const_cast<char*>(FSpath_string.c_str());
    face_coordinates_vector.clear();
    face_coordinates_vector.resize(CellNumbs.at(2), make_tuple(0,0,0));
//REPAIR    cout << "face vector size: " << face_coordinates_vector.size() << endl;

    if (is_file_exists(FSpath_string)) {
        face_coordinates_vector = TuplesReader(FSpath);
    }
    else {
        cout << "Calculation of face coordinates started now by find_aGBseed function" << endl;
        for (unsigned int fnumber = 0; fnumber < CellNumbs.at(2)-1; ++fnumber) {
            face_coordinates_vector.at(fnumber) = find_aGBseed(fnumber, paths, CellNumbs, grain_coordinates);
//REPAIR
     if(fnumber % 10 == 0)  cout << "Total # of faces " << CellNumbs.at(2)<< " face # " << fnumber << " "<< get<0>(face_coordinates_vector.at(fnumber)) << " "
                 << get<1>(face_coordinates_vector.at(fnumber)) << " " << get<2>(face_coordinates_vector.at(fnumber)) << endl;
        }
        cout << "Calculations was done!" << endl;
        ofstream OutFaceCoord; OutFaceCoord.open(FSpath, ios::trunc);
        if(OutFaceCoord.is_open())
            for (auto fcoord : face_coordinates_vector)
                OutFaceCoord << get<0>(fcoord) << " " << get<1>(fcoord) << " " << get<2>(fcoord) << endl;
    }

    /// All vertex coordinates
    string VCpath_string = input_folder + "vertex_seeds.txt"s;
    char* VCpath = const_cast<char*>(VCpath_string.c_str());
    vertex_coordinates_vector.resize(CellNumbs.at(0), make_tuple(0,0,0));
    vertex_coordinates_vector = TuplesReader(VCpath);
//REPAIR    cout << "vertex_coordinates_vector size: " << vertex_coordinates_vector.size() << endl; //exit(0); exit(0);

    /// All grain volumes
    string GVpath_string = input_folder + "grain_volumes.txt"s;
    char* GVpath = const_cast<char*>(GVpath_string.c_str());
    grain_volumes_vector.resize(CellNumbs.at(3));
    grain_volumes_vector = dVectorReader(GVpath);

    /// All face areas
    string GBApath_string = input_folder + "face_areas.txt"s;
    char* GBApath = const_cast<char*>(GBApath_string.c_str());
    face_areas_vector.resize(CellNumbs.at(3));
    face_areas_vector = dVectorReader(GBApath);

/// Special feature for the 2D case ("grains" -> 2-cells; "faces" -> 1-cells; "triple junctions (lines) -> 0-cells" in 2D):
// Very important!
// It needs to make everything similar in the code for 2D and 3D cases with the same output of the DCC Generator tool (https://github.com/PRISBteam/Voronoi_DCC_Analyser) software (the same set of A0, A1..., B2 matrices)
    if (dim == 2) { CellNumbs.push_back(CellNumbs.at(2)); CellNumbs.at(2) = CellNumbs.at(1); CellNumbs.at(1) = CellNumbs.at(0); CellNumbs.at(0) = NULL; }

/// CellNumbs output
    cout << "=====================================================================================" << endl;
    unsigned int t_length = 0;
    for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl;
    //---------------------------------------------------------------------------
/// The PRIMAL DATA STRUCTURES and concepts - special (SFS) and ordinary (OFS) face sequences, defect state vectors (DSVs) and defect configurations (DCs)
    // Order of newly generated special faces | Process time
    std::vector <unsigned int> special_faces_sequence, ordinary_faces_sequence, current_sfaces_sequence, current_ofaces_sequence, crack_faces_sequence, current_cracks_sequence; // Variable sequences (in order of their generation) of special and ordinary Faces and Cracks
    vector<unsigned int> sub_sfaces_sequence, frac_sfaces_sequence;
    // State vector | Special faces IDs
    vector <unsigned int> State_sVector(CellNumbs.at(2), 0), current_State_sVector(CellNumbs.at(2), 0), State_cVector(CellNumbs.at(2), 0), current_State_cVector(CellNumbs.at(2), 0);
    // MAX fraction of special Faces | Calculation limit
    max_sFaces_fraction = ConfigVector.at(2), max_cFaces_fraction = ConfigVector.at(3);

    cout << "=========================================================================" << endl << "\t\t\t\t\t[\tStart of the DCC Processing Tool\t]\t\t\t\t\t" << endl << "-------------------------------------------------------------------------" << endl;

/// ==========================================================================================================================================
/// ================================================= THE LIST MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
/// In the LIST mode all the functions are calling one after another without loops
    double P_time = 0, K_time = 0, C_time = 0, W_time = 0;
    double mu_f = 0, sigm_f = 0; // only for L type simulations (!)
    vector<vector<int>> RW_series_vector; // only for L type simulations (!)

/// Clearance of files
//REPAIR    OfStreams_trancator();
    special_face_design.clear(); // clearing file

    /// THE CURRENT TASKS ARE PLACED HERE ///
/// ==========================================================================================================================================

//#include "tasks/task_macrocrack.cpp"

/// ==========================================================================================================================================

/// I: DCC_Processing module
//    if (ProcessingON(confpath, time_step_one)) { // if DCC_Processing is SWITCH ON in the config.txt file
    cout << "START of the DCC Processing module" << endl;
/*        if (P_type != "L") {
            special_face_design = DCC_Processing(special_faces_sequence, State_sVector, P_type);
        } else { // Vector of sfaces strips generated by RW
            for (int mu = 5; mu < 6; mu++) { //mu_f_max = 20
                mu_f = (double) mu;
                for (int sigm = 3; sigm < 4 ; sigm++) { //sigm_f_max = 0.7
                    sigm_f = (1.0 / (double) sigm) * 0.9;
                    special_face_design = DCC_Processing(special_faces_sequence, State_sVector, P_type, mu_f, sigm_f, RW_series_vector);
                } //  for (int sigm = 1; sigm < 10 ; sigm++)
            } // for (int mu = 2; mu < 3; mu++)
        }*/

    special_face_design = DCC_Processing(special_faces_sequence, State_sVector, "R");

/// ===== Elapsing time Processing ================
//        if (ProcessingON(confpath, time_step_one)) {
    unsigned int Processing_time = clock();
    P_time = (double) Processing_time;
    cout << "Processing time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl;
    cout << "-------------------------------------------------------" << endl;
//       } // end if(ProcessingON(confpath))

/// Cut special face sequence
    double s_fraction = 0.1; // variable parameter
    frac_sfaces_sequence.clear(); // part of sequence with certain fraction
    for (unsigned int utr = 0; utr < special_faces_sequence.size()*s_fraction; ++utr) {
        frac_sfaces_sequence.push_back(special_faces_sequence.at(utr));
    }
    cout << "frac_sfaces_sequence size:  " << frac_sfaces_sequence.size() << endl;

/// 0: DCC_Subcomplex module
//    if (SubcomplexON(confpath, time_step_one)) { // if DCC_Processing is SWITCH ON in the config.txt file
    cout << "START of the DCC Subcomplex module" << endl;
    unsigned int subcomplex_id_new = 0;
    subcomplex new_Cut(subcomplex_id_new);
//       if (S_type == "P") { // subcomplex type  "plane cut"
    new_Cut = subcomplex(subcomplex_id_new);

  //  int D_initial_number = 1, nDsteps = 100;
 //   for (int Dp_number = D_initial_number; Dp_number < nDsteps; ++ Dp_number) {
 //       double Dp_value = (1.0/(double) nDsteps) * (double) Dp_number;
 //       double a_p = 0.0, b_p = 0.0, c_p = 1.0, D_p = Dp_value;
 //       DCC_Subcomplex(new_Cut, frac_sfaces_sequence, sub_sfaces_sequence, crack_faces_sequence, a_p, b_p, c_p, D_p);

        double a_p = 0.0, b_p = 0.0, c_p = 1.0, D_p = 0.5; /// Must be the same as in the DCC_Subcomplex!!!
    DCC_Subcomplex(new_Cut, frac_sfaces_sequence, sub_sfaces_sequence, crack_faces_sequence, a_p, b_p, c_p, D_p);
//            }
//    } // end of if (SubcomplexON(confpath, time_step_one))

    cout << "plane_cut by DCC_Subcomplex was successfully created!" << endl;
/// ===== Elapsing time Subcomplex ================
    unsigned int Subcomplex_time = clock();
    double Sub_time = (double) Subcomplex_time - P_time;
    cout << "Subcomplex time is equal to  " << Sub_time / pow(10.0, 6.0) << "  seconds" << endl;
    cout << "-------------------------------------------------------" << endl;

    ///////////////// LOOP over crack steps /////////////////////
    OutCrackEnergies_file.open(output_folder + "Crack_energies.txt"s, ios::trunc); // ofstream for crack energies
    OutCrackEnergies_file << "(1) sfaces fraction" << "  " << "(2) Crack length ratio" << "  " << "(3) external vonMizes stress" << "  " << "(4) macro_crack.surface_energy" << "  " << "(5) macro_crack.bridging_energy" << "  " << "(6) multiple_cracking_energy" << "  " << "(7) crack_fraction" << endl;

    int crack_id_new = 0;
    int nsteps = 10;
    int initial_step = 1;
    for (int crack_length_number = initial_step; crack_length_number < nsteps; ++crack_length_number) {

/// Half-plane subcomplex
    subcomplex half_plane_cut(new_Cut);
    half_plane_cut.a_n = a_p; half_plane_cut.b_n = b_p; half_plane_cut.c_n = c_p; half_plane_cut.D_plane = D_p;
    double length_ratio = (1.0/(double) nsteps) * (double) crack_length_number;
    double crack_length = length_ratio; //* Lx_size; // Lx_size here Lx only now!!!
    int half_subcomplex_id_new = 0;
// half_plane_cut

    /// Half plane
    half_plane_cut = new_Cut.Get_half_plane(new_Cut, crack_length, sub_sfaces_sequence); // half-plane method for a subcomplex

///News
    cout << "half_plane_cut by Get_half_plane(plane_cut, crack_length) \"subcomplex\" class method was successfully created!" << endl;
    unsigned int HalfPlane_time = clock();
    double HP_time = (double) HalfPlane_time - Sub_time;
    cout << "half plane cut takes time: " << HP_time/ pow(10.0, 6.0) << " s" << endl;

/// Multiphysics module : Macrocrack and the vector of macrocracks

/// number of "steps" in a crack growth path
    ++crack_id_new; // new crack id
    // vector of macro-cracks as a PCC objects
    large_cracks_vector.push_back(macrocrack(crack_id_new, half_plane_cut));
    large_cracks_vector.at(crack_length_number-initial_step).Set_crack_plane(); // subcomplex based

    ///current crack
    macrocrack large_crack = large_cracks_vector.at(crack_length_number-initial_step);

    /// External stress value
    double external_vonMizes_stress = 3.0 * pow(10.0, 8.0); // 300 MPa
    face_elastic_energies.resize(CellNumbs.at(2), 0.0);

    face_elastic_energies = DCC_Multiphysics(large_crack, external_vonMizes_stress);

    ///News
    cout << "face_elastic_energies vector by DCC_Multiphysics module was successfully calculated!" << endl;
    unsigned int MultiPhys = clock();
    double MP_time = (double) MultiPhys - HP_time;
    cout << "DCC_Multiphysics module takes time: " << MP_time / pow(10.0, 6.0) << " s" << endl;

/// III: DCC_Kinetic module
//if (KineticON(confpath, time_step_one)) { // if DCC_Kinetic is SWITCH ON in the config.txt file

    cout << "START of the DCC Kinetic module" << endl;
    kface_sequence = DCC_Kinetic(face_elastic_energies, frac_sfaces_sequence, large_crack, "F"s);

/// ===== Elapsing time Kinetic ================
    unsigned int Kinetic_time = clock();
    K_time = (double) Kinetic_time - MP_time;
    cout << "Kinetic time is equal to  " << K_time / pow(10.0, 6.0) << "  seconds" << endl;
    cout << "-------------------------------------------------------" << endl;
//}// end if(KineticON(confpath))

    cout << "START of the DCC Writer module" << endl;
// DCC_Writer(special_faces_sequence, kface_sequence, State_sVector, mu_f, sigm_f, RW_series_vector, P_type);
//DCC_subcomplex_Writer(subcomplex_id_new, plane_cut, half_subcomplex_id_new, half_plane_cut, special_faces_sequence);
    DCC_subcomplex_Writer(0, new_Cut, 1, half_plane_cut, frac_sfaces_sequence, kface_sequence, large_crack,
                          external_vonMizes_stress);

/// ===== Elapsing time Writer ================
    unsigned int Writer_time = clock();
    W_time = (double) Writer_time - K_time;
//    W_time = (double) Writer_time - K_time - MP_time - HP_time - Sub_time - P_time;
    cout << "------" << endl;
    cout << "Writer time is equal to  " << W_time / pow(10.0, 6.0) << "  seconds" << endl;

} // end of for loop over crack lengths
// } // end of for loop over the planes Dp
    OutCrackEnergies_file.close();


    /*
    if ( SIMULATION_MODE(confpath) == "SIMULATION MODE(LIST)"s) { /// SIMULATION MODE: #LIST
// Clear S_Vector //std::fill(S_Vector.begin(), S_Vector.end(), 0);

/// I: DCC_Processing module
    if (ProcessingON(confpath, time_step_one)) { // if DCC_Processing is SWITCH ON in the config.txt file
        cout << "START of the DCC Processing module" << endl;
        if (P_type != "L") {
            special_face_design = DCC_Processing(special_faces_sequence, State_sVector, P_type);
            } else { // Vector of sfaces strips generated by RW
            for (int mu = 5; mu < 6; mu++) { //mu_f_max = 20
                mu_f = (double) mu;
                for (int sigm = 3; sigm < 4 ; sigm++) { //sigm_f_max = 0.7
                    sigm_f = (1.0 / (double) sigm) * 0.9;
                    special_face_design = DCC_Processing(special_faces_sequence, State_sVector, P_type, mu_f, sigm_f, RW_series_vector);
                } //  for (int sigm = 1; sigm < 10 ; sigm++)
            } // for (int mu = 2; mu < 3; mu++)
        }

/// ===== Elapsing time Processing ================
        if (ProcessingON(confpath, time_step_one)) {
            unsigned int Processing_time = clock();
            P_time = (double) Processing_time;
            cout << "Processing time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl;
            cout << "-------------------------------------------------------" << endl;
        } // end if(ProcessingON(confpath))

        /// Loop over special_faces_sequence with the step = Face_fraction/ number_of_fsteps
        long number_of_fsteps = 100, number_of_csteps = 1; /// NUMBER OF STEPS for loops over the fractions of special Faces and Cracks
        double Face_fraction = 0.01, Crack_fraction = 0.00; /// INITIAL fractions of special Faces and Cracks
        double dp = 0, df = 0; //increments (dp - special Faces fraction, df - cracked faces fraction)

        /// Writer on the screen
        cout << "************** ======== LOOPS OVER SPECIAL AND CRACKE FACES ======== **************"s << endl;
        cout << "Number of steps over Special Faces loop = " << number_of_fsteps << endl;
        cout << "Number of steps over Cracked Faces loop = " << number_of_csteps << endl;
        cout << "Initial fraction of Special Faces in the loop = " << Face_fraction << endl;
        cout << "Initial fraction of Cracked Faces in the loop = " << Crack_fraction << endl;
        cout << "Maximal fraction of Special Faces in the loop = " << max_sFaces_fraction
             << endl; // is taken from config.txt file
        cout << "Maximal fraction of Cracked Faces in the loop = " << max_cFaces_fraction
             << endl; // is taken from config.txt file
        cout << "==================================================================================="s << endl;

        /// ofstreams for Face Conductivity and Face Laplacian output are opening again with "app" mode before the loops
        // all these ofstreams will be sent to the Characterisation function for the specific output
///                        OutElCondfile.open(FCond_dir, ios::app);

        /// ======== The first loop over all Face_fraction with the step dp ===========>>>
        if (max_sFaces_fraction > 0) {
            for (long i = 0; i < number_of_fsteps; ++i) { // loop over all Face_fraction
                Face_fraction += dp;
                unsigned int special_Faces_numb = Face_fraction * CellNumbs.at(2);

                /// creation of the current_sfaces_sequence as the part (before special_Faces_numb) of the special_faces_sequence
                unsigned int itu = 0;
                current_sfaces_sequence.clear();
                for (unsigned int sfe: special_faces_sequence) {
                    if (itu <= special_Faces_numb) current_sfaces_sequence.push_back(sfe);
                    ++itu;
                } // REPAIR for (auto sfe : current_sfaces_sequence) cout << sfe << "\t"; cout << endl;

                /// creation of the current_State_Vector by current_sfaces_sequence
                unsigned int itr = 0;
                fill(current_State_sVector.begin(), current_State_sVector.end(), 0);
                for (auto sv: current_sfaces_sequence) current_State_sVector.at(sv) = 1;

/// II_: DCC_Multiphysics module and others

                /// plane subcomplex
                unsigned int subcomplex_id_new = 0;
                subcomplex plane_cut(subcomplex_id_new);
                DCC_Subcomplex(plane_cut, special_faces_sequence, crack_faces_sequence); // output subcomplex type plane_cut

                /// News
                cout << "plane_cut by DCC_Subcomplex was successfully created!" << endl;
                unsigned int s1_time = clock();
                s1_time = (double) s1_time - P_time;
                cout << "plane_cut takes time: " << s1_time/ pow(10.0, 6.0) << " s" << endl;

                /// half-plane subcomplex
                double length_ratio = 0.1;
                double crack_length = length_ratio * Lx_size; // Lx_size here Lx only now!!!
                int half_subcomplex_id_new = 0;
                subcomplex half_plane_cut(half_subcomplex_id_new);
                half_plane_cut = plane_cut.Get_half_plane(plane_cut, crack_length); // half-plane method for a subcomplex

                ///News
                cout << "half_plane_cut by Get_half_plane(plane_cut, crack_length) \"subcomplex\" class method was successfully created!" << endl;
                unsigned int s2_time = clock();
                s2_time = (double) s2_time - s1_time;
                cout << "plane_cut takes time: " << s2_time/ pow(10.0, 6.0) << " s" << endl;

                DCC_subcomplex_Writer(subcomplex_id_new, plane_cut, half_subcomplex_id_new, half_plane_cut, special_faces_sequence);

                /// Macrocrack and the vector of macrocracks
                int crack_id_new = 0;
                // vector of macro-cracs as a PCC objects
                large_cracks_vector.push_back(macrocrack(crack_id_new, plane_cut));

                /// External stress value
                double external_vonMizes_stress = pow(10.0,8.0); // 100 MPa
                face_elastic_energies = DCC_Multiphysics(large_cracks_vector, external_vonMizes_stress);

                ///News
                cout << "face_elastic_energies vector by DCC_Multiphysics module was successfully calculated!" << endl;
                unsigned int s3_time = clock();
                s3_time = (double) s3_time - s2_time;
                cout << "plane_cut takes time: " << s3_time/ pow(10.0, 6.0) << " s" << endl;

/// III: DCC_Kinetic module
                if (KineticON(confpath, time_step_one)) { // if DCC_Kinetic is SWITCH ON in the config.txt file

                    cout << "START of the DCC Kinetic module" << endl;
                    kface_sequence = DCC_Kinetic(face_elastic_energies, special_faces_sequence, K_type);

/// ===== Elapsing time Kinetic ================
                    unsigned int Kinetic_time = clock();
                    K_time = (double) Kinetic_time - P_time;
                    cout << "Kinetic time is equal to  " << K_time / pow(10.0, 6.0) << "  seconds" << endl;
                    cout << "-------------------------------------------------------" << endl;
                }// end if(KineticON(confpath))

/// II: DCC_Characterisation module
                if (CharacterisationON(confpath, time_step_one)) {

                    cout << "--------------------------------" << endl;
                    cout << "START of the DCC Structure Characterisation module" << endl;
                    DCC_StructureCharacterisation(State_sVector, special_faces_sequence, crack_faces_sequence,
                                                  ConfigVector);

/// ===== Elapsing time Writer ================
                    unsigned int Characterisation_time = clock();
                    C_time = (double) Characterisation_time - P_time - K_time;

                    cout << "Characterisation time is equal to  " << C_time / pow(10.0, 6.0) << "  seconds"
                         << endl;
                    cout << "-------------------------------------------------------" << endl;
                }// end if(CharacterisationON(confpath))

                if (WriterON(confpath, time_step_one)) { // simply output ON/OFF for the DCC_Writer module on the screen
                    cout << "START of the DCC Writer module" << endl;
                    DCC_Writer(special_faces_sequence, kface_sequence, State_sVector, mu_f, sigm_f, RW_series_vector, P_type);

/// ===== Elapsing time Writer ================
                    unsigned int Writer_time = clock();
                    double W_time = (double) Writer_time - P_time - K_time - C_time;

                    cout << "------" << endl;
                    cout << "Writer time is equal to  " << W_time / pow(10.0, 6.0) << "  seconds" << endl;
                } // end if(WriterON(confpath))

///                    } // end of for (int mu = 1; mu < 20; mu++)
///                } // end of for (int sigm = 1; sigm < 20; sigm++)
                //                   }  // end of else -> if(P_type != "L")

                /// Increment of special faces fraction
                dp = max_sFaces_fraction / (double) number_of_fsteps; // Face fraction increment
            } /// END of loop : for (number_of_fsteps)
        } //if max_sFaces_fraction > 0
///        } //ProcessingON
    else cout << "WARNING!!! max_sFaces_fraction = 0" << endl;

} /// end SIMULATION MODE else (LIST)
*/

/// ===== Elapsing time ================>
    unsigned int end_time = clock();
    double fulltime = (double) end_time;
    cout << "-------------------------------------------" << endl;
    cout << dim << "D " << "runtime is equal to  " << fulltime/pow(10.0,6.0) <<  "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the VoroCAnalyser program\t]\t\t\t\t\t" << endl << "=========================================================================" << endl;

    return 0;
} /// The END of Main function

/// ================================== DEFINED IN MAIN FUNCTIONS ==================================///

/// ================== Initial configuration - reading and output =================>
string SIMULATION_MODE(char* config) {
    std::string line;
    ifstream inConf(config);
    string typeCharacterisationON;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
// REPAIR           cout << line << endl;
            if (!line.compare("SIMULATION MODE(FRACTURE_LOOP)"s))
            { typeCharacterisationON = "SIMULATION MODE(FRACTURE_LOOP)"s; cout << "FRACTURE_LOOP SIMULATION MODE"s << endl; }
            else if (!line.compare("SIMULATION MODE(INDEX_LOOP)"s))
            { typeCharacterisationON = "SIMULATION MODE(INDEX_LOOP)"s; cout << "INDEX_LOOP SIMULATION MODE"s << endl; }
            else if (!line.compare("SIMULATION MODE(ISING_LOOP)"s))
            { typeCharacterisationON = "SIMULATION MODE(ISING_LOOP)"s; cout << "ISUNG_LOOP SIMULATION MODE"s << endl; }
            else if (!line.compare("SIMULATION MODE(LIST)"s))
            { typeCharacterisationON = "SIMULATION MODE(LIST)"s; cout << "LIST SIMULATION MODE"s << endl; }
    } else cout << "SIMULATION_MODE() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    return typeCharacterisationON;
} // end string SIMULATION_MODE(char* config)

bool ProcessingON(char* config, bool time_step_one) {
    std::string line;
    ifstream inConf(config);
    bool isProcessingON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
// REPAIR            cout << line << endl;
            if (!line.compare("DCC_Processing SWITCHED ON"s)) {
                isProcessingON = 1;
                if (time_step_one == 1) cout << "ON    | DCC_Processing"s << endl;
                return isProcessingON;
            }
    } else cout << "ProcessingON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "OFF    | DCC_Processing"s << endl;
    return isProcessingON;
}

bool KineticON(char* config, bool time_step_one) {
    std::string line;
    ifstream inConf(config);
    bool isKineticON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
// REPAIR            cout << line << endl;
            if (!line.compare("DCC_Kinetic SWITCHED ON"s)) {
                isKineticON = 1;
                if (time_step_one == 1) cout << "ON   | DCC_Kinetic"s << endl;
                return isKineticON;
            }
    } else cout << "KineticON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "OFF   | DCC_Kinetic"s << endl;
    return isKineticON;
}

bool CharacterisationON(char* config, bool time_step_one) {
    std::string line;
    ifstream inConf(config);
    bool isCharacterisationON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
// REPAIR            cout << line << endl;
            if (!line.compare("DCC_Characterisation SWITCHED ON"s)) {
                isCharacterisationON = 1;
                if (time_step_one == 1)    cout << "ON   | DCC_Characterisation"s << endl;
                return isCharacterisationON;
            }
    } else cout << "isCharacterisationON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "OFF   | DCC_Characterisation"s << endl;
    return isCharacterisationON;
}

bool WriterON(char* config, bool time_step_one) {
    std::string line;
    ifstream inConf(config);
    bool isWriterON = 0;

    if (inConf) { // If the file was successfully open, then
        while(getline(inConf, line, '\n'))
// REPAIR            cout << line << endl;
            if (!line.compare("DCC_Writer SWITCHED ON"s)) {
                isWriterON = 1;
                if (time_step_one == 1) cout << "ON    | DCC_Writer"s << endl;
                return isWriterON;
            }
    } else cout << "WriterON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "OFF   | DCC_Writer"s << endl;
    return isWriterON;
}

    bool SubcomplexON(char* config, bool time_step_one) {
        std::string line;
        ifstream inConf(config);
        bool isSectionON = 0;

        if (inConf) { // If the file was successfully open, then
            while(getline(inConf, line, '\n'))
// REPAIR            cout << line << endl;
                if (line.compare("DCC_Section SWITCHED ON"s)) {
                    isSectionON = 1;
                    if (time_step_one == 1) cout << "ON    | DCC_Section"s << endl;
                    return isSectionON;
                }
        } else cout << "SubcomplexON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

        if (time_step_one == 1) cout << "OFF   | DCC_Section"s << endl;
        return isSectionON;
    }

/// Reading and Output of the configuration file
std::vector<double> confCount(char* config, string &Subcomplex_type, string &Processing_type, string &Kinetic_type, string &input_folder, string &output_folder) {
    vector<double> res;

    string line;
    double p_max = 0, f_max = 0;

    ifstream inConf(config);
    if (inConf) { // If the file was successfully open, then
        while (getline(inConf, line, '\n')) {
            // Number of types
            for (auto it: line) {
                if (it == '}') Subcomplex_type = line.at(0); // simulation Section type
                else if (it == '&') Processing_type = line.at(0); // simulation Processing type
                else if (it == '`') Kinetic_type = line.at(0); // simulation Kinetic type

                else if (it == '@') res.push_back(line.at(0) - '0'); // dimension of the problem; res[1]
                else if (it == '!') res.push_back(line.at(0) - '0'); // number of i(X) designs; res[0]
                    ///?? does it properly working with 10, 100 etc
                    //else if (it == '!') res.push_back(line.at(0) - '0'); // number of special Face types; res[0]

                else if (it == '?') {
                    stringstream line1_stream(line);
                    line1_stream >> p_max;
                    res.push_back(p_max);
                } // MAX fraction of Faces (calculation limit); res[2]
                else if (it == '^') {
                    stringstream line2_stream(line);
                    line2_stream >> f_max;
                    res.push_back(f_max);
                } // MAX fraction of Cracks (calculation limit); res[3]
                else if (it == '~') {
                    stringstream line3_stream(line);
                    line3_stream >> input_folder;
                } // input folder path input_folder = const_cast<char*>(input.c_str());
                else if (it == '$') {
                    stringstream line4_stream(line);
                    line4_stream >> output_folder;
                } // output folder path input_folder = const_cast<char*>(input.c_str());

                else if (it == '#') res.push_back(1); // 1 and # means accept - the parameter will be calculated
                else if (it == '%') res.push_back(0); // 0 and % means ignore - the parameter will not be calculated
            }
        }

    } else cout << "The file " << config << " cannot be read" << endl; // If something goes wrong

    res.size();
    cout << "The problem dimension that is the maximum dimension k_max of k-Cells\t | "s << res.at(0) << endl;
    cout << "The number of i(X) designs (0 - random, 1 - random + smax,.. )\t\t\t | " << res.at(1) << endl;
    if (SubcomplexON(config, 0)) cout << "Subcomplex type ('P' or 'H')\t\t\t\t\t\t\t\t\t | "s << Subcomplex_type << endl;
    if (ProcessingON(config, 0)) cout << "Processing type ('R', 'S', 'C' or 'X')\t\t\t\t\t\t\t\t\t | "s << Processing_type << endl;
    //cout << "Calculation type ('R', 'RW', 'S', 'F', 'I' or 'X'):\t\t\t\t\t\t | "s << Processing_type << endl;
    //cout << "The number of special Face (2-cells) types:\t\t\t\t\t\t\t\t | " << res.at(1) << endl;
    cout << "MAX fraction of Faces (calculation limit) \t\t\t\t\t\t\t\t | " << res.at(2) << endl;
    if (KineticON(config, 0)) cout << "Kinetic type ('W' or 'I')\t\t\t\t\t\t\t\t\t\t\t\t | "s << Kinetic_type << endl;
    if (KineticON(config, 0)) cout << "MAX fraction of Cracks (calculation limit)\t\t\t\t\t\t\t\t | " << res.at(3) << endl;
    cout << endl;
    cout << "Input folder:\t" << input_folder << endl;
    cout << "Output folder:\t" << output_folder << endl;
    cout << endl;
    if (CharacterisationON(config, 0)) {
        cout << "Nodes types statistics, indices and configuration entropy:     "; if (res.at(4) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Edges types statistics, indices and configuration entropy:     "; if (res.at(5) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Faces types statistics and structural indices:                 "; if (res.at(6) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Grain types statistics, indices and configuration entropy:     "; if (res.at(7) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Node Laplacian with its spectrum for the Nodes network:       "; if (res.at(8) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Edge Laplacian with its spectrum for the Edges network:       "; if (res.at(9) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Face  Laplacian with its spectrum for the Faces network:       "; if (res.at(10) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Grain Laplacian with its spectrum for the Grains network:      "; if (res.at(11) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Tutte polynomial for the special network:                      "; if (res.at(12) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    }

    return res;
}

//Erase First Occurrence of given  substring from main string
void eraseSubStr(std::string & mainStr, const std::string & toErase)
{
    // Search for the substring in string
    size_t pos = mainStr.find(toErase);
    if (pos != std::string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

/// Archive
/* !!!!!!!! WORKING CODE !!!!!!!!!
 * /// ===========================================================================================================================================
/// ============================================= THE FRACTURE LOOP MODE STARTS HERE ===================================================================
/// ===========================================================================================================================================
    else if ( SIMULATION_MODE(confpath) == "SIMULATION MODE(FRACTURE_LOOP)"s) { /// SIMULATION MODE: #LOOP
        ///  ====== Analytical solutions ========>
        // The function from TJsLab.h giving the analytical solution in Random and quasi-random (with crystallographical restrictions) cases for HAGBs fractions
       // TJsAnalytics(1000, odir);
//         TJsAnalytics_indices(2.0, 1.5, 1.7); exit(51);

            /// ====== Processing ========>
        if (ProcessingON(confpath, time_step_one)) {
            cout << "START of DCC Processing Module" << endl;
//            DCC_Processing3D(State_Vector, special_faces_sequence, Processing_type, max_sFaces_fraction, number_of_types, CellNumbs, paths, odir, ip_index); // a function from Processing_Functions.h
              DCC_Processing3D( State_sVector, special_faces_sequence, P_type);

        cout << "The DCC_Processing successfully finished with the output of special_faces_sequence and State_Vector"s << endl;
        cout << endl;
        }
        else {
            /// ====== Reading from file instead of DCC_Processing3D ( ) if it is off ============>
         string SFS_path = output_folder + "Smax/SpecialGrainBoundaries.txt"s;
         char* SFS_dir = const_cast<char *>(SFS_path.c_str());
            special_faces_sequence = Smax_sequence_reader(SFS_dir); // a function from Processing_Functions.h
            for (auto itr: special_faces_sequence) State_sVector.at(itr) = 1;
        } // end of else

        /// Loop over special_faces_sequence with the step = Face_fraction/ number_of_fsteps
    long number_of_fsteps = 30, number_of_csteps = 30; /// NUMBER OF STEPS for loops over the fractions of special Faces and Cracks
    double Face_fraction = 0.01, Crack_fraction = 0.01; /// INITIAL fractions of special Faces and Cracks
    double dp = 0, df = 0; //increments (dp - special Faces fraction, df - cracked faces fraction)

   /// Writer on the screen
    cout << "************** ======== LOOPS OVER SPECIAL AND CRACKE FACES ======== **************"s << endl;
    cout << "Number of steps over Special Faces loop = " << number_of_fsteps << endl;
    cout << "Number of steps over Cracked Faces loop = " << number_of_csteps << endl;
    cout << "Initial fraction of Special Faces in the loop = " << Face_fraction << endl;
    cout << "Initial fraction of Cracked Faces in the loop = " << Crack_fraction << endl;
    cout << "Maximal fraction of Special Faces in the loop = " << max_sFaces_fraction << endl; // is taken from config.txt file
    cout << "Maximal fraction of Cracked Faces in the loop = " << max_cFaces_fraction << endl; // is taken from config.txt file
    cout << "==================================================================================="s << endl;

        /// Output (ofstream) paths for different files from Characterisation module
        string output_TJs_dir = output_folder + "TJsLab_TJsTypes.txt"s; char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str());  // From string to char for the passing folder path to a function
        string output_Entropy_dir = output_folder + "TJsLab_TJsEntropy.txt"s; char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str());
        string Face_Cond_dir = output_folder + "Face_Conductivity.txt"s; char* FCond_dir = const_cast<char*>(Face_Cond_dir.c_str());
        string Face_Laplacian_dir = output_folder + "SpecialFacesLaplacian.txt"s; char* FLap_dir = const_cast<char*>(Face_Laplacian_dir.c_str());

        /// Creation of the two files with TJs fraction and Configuration entropies output
        // ofstream for entropies OutSFile
        ofstream OutSFile; OutSFile.open(cEntropy_dir, ios::trunc);
        OutSFile << "\t(1)Fraction_of_special_Faces\t" << " " << "\t(2)Fraction_of_cracks\t" << " " << "\t(3)Von-Neumann_Special_Faces_entropy\t" << " " << "\tj(4)0_TJs\t"<< " " << "\t(5)j1_TJs\t"<< " " << "\t(6)j2_TJs\t"<< " " << "\t(7)j3_TJs\t" << " " << "\t(8)p*=j1+2*j2+3*j3/3\t" << " " << "\t(9)Special_Faces_Entropy\t" << " " << "\t(10)Special_Faces_median_entropy\t" << " " << "\t(11)Special_Faces_skrew_entropy\t" << " " << "\t(12)Special_Faces_informativeness\t" << " "
                    << "\t(13)Von-Neumann_Cracked_Faces_entropy\t" << " " << "\t(14)jc0_cTJs\t"<< " " << "\t(15)jc1_cTJs\t"<< " " << "\t(16)jc2_cTJs\t"<< " " << "\t(17)jc3_cTJs\t" << " " << "\t(18)f*=jc1+2*jc2+3*jc3/3\t"<< " " << "\t(19)Special_Cracks_Entropy\t" << " " << "\t(20)Cracked_Face_median_entropy\t" << " " << "\t(21)Cracked_Face_skrew_entropy\t" << " " << "\t(22)Cracked_Face_informativeness\t" << endl;
        OutSFile << endl;
        OutSFile.close();
        // ofstream for TJs OutTJsFile
        ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::trunc);
        OutTJsFile.close();

        /// Output to file Special Face Laplacian eigenvalues
        // Special Face Conductivity output
        ofstream OutElCondfile; OutElCondfile.open(FCond_dir, ios::trunc);
        OutElCondfile << "\t(1)Fraction_of_special_Faces\t" << " " << "\t(2)Fraction_of_cracks\t" << " " << "\t(3)Direct_sum_of_Special_Face_Laplacian_eigenvalues\t"<< " " << "\t(4)Special_Face_conductivity\t" << " " << "\t(5)Number_of_special_Face_components\t" << " " <<"\t(6)Special Face informativeness\t" << " "
                      <<"\t(7)Direct_sum_of_Cracked_Face_Laplacian_eigenvalues\t" << " " << "\t(8)Cracked_Face_conductivity\t" << " " << "\t(9)Number_of_cracked_Face_components\t" << " " << "\t(10)Cracked Face informativeness\t" << endl;
        OutElCondfile << endl;
        OutElCondfile.close();
        // Special Face Laplacian output
        ofstream OutFLfile; OutFLfile.open(FLap_dir, ios::trunc);
        OutFLfile.close();

        /// ofstreams for Face Conductivity and Face Laplacian output are opening again with "app" mode before the loops
        // all these ofstreams will be sent to the Characterisation function for the specific output
        OutElCondfile.open(FCond_dir, ios::app);
        OutFLfile.open(FLap_dir, ios::app);

        /// ======== The first loop over all Face_fraction with the step dp ===========>>>
    if (max_sFaces_fraction > 0) {
        for (long i = 0; i < number_of_fsteps; ++i) { // loop over all Face_fraction
            Face_fraction += dp;
            unsigned int special_Faces_numb = Face_fraction * CellNumbs.at(2);

            /// creation of the current_sfaces_sequence as the part (before special_Faces_numb) of the special_faces_sequence
            unsigned int itu = 0;
            current_sfaces_sequence.clear();
            for (unsigned int sfe: special_faces_sequence) {
                if (itu <= special_Faces_numb) current_sfaces_sequence.push_back(sfe);
                ++itu;
            } // REPAIR for (auto sfe : current_sfaces_sequence) cout << sfe << "\t"; cout << endl;

            /// creation of the current_State_Vector by current_sfaces_sequence
            unsigned int itr = 0;
            fill(current_State_sVector.begin(), current_State_sVector.end(), 0);
            for (auto sv: current_sfaces_sequence) current_State_sVector.at(sv) = 1;

            /// ====== Fracture Kinetic module ========>
            if (KineticON(confpath, time_step_one) && max_cFaces_fraction > 0)
//                crack_faces_sequence = DCC_Kinetic(current_sfaces_sequence, K_type);
// REPAIR        for (auto cn : crack_sequence)  cout << cn << "\t"; cout << endl; //        exit(19);

            /// ================== The second loop over all crack_Faces_numb ===================================>>>
        if (max_cFaces_fraction > 0) {
            for (long y = 0; y < number_of_csteps; ++y) { /// loop over all Cracks
                Crack_fraction += df;

                unsigned int crack_Faces_numb = Crack_fraction * CellNumbs.at(2);
                unsigned int itc = 0;
                current_cracks_sequence.clear();
                for (unsigned int sfe: crack_faces_sequence) {
                    if (itc <= crack_Faces_numb) current_cracks_sequence.push_back(sfe);
                    ++itc;
                } // REPAIR       for (auto sfe : current_cracks_sequence) cout << sfe << "\t"; cout << endl;

            /// ====== Structure Characterisation module ========>
            if (CharacterisationON(confpath, time_step_one)) {
//////////               DCC_StructureCharacterisation(current_State_sVector, current_sfaces_sequence, current_cracks_sequence, ConfigVector, CellNumbs, paths, odir, OutFLfile, OutElCondfile, OutSFile);
//                 DCC_StructureCharacterisation(current_State_sVector, current_sfaces_sequence, current_cracks_sequence, ConfigVector);

                } else break;

                // the following time steps
                time_step_one = 0;
                /// Increment of crack fraction
                df = max_cFaces_fraction / (double) number_of_csteps;
            } /// END of loop : for (number_of_csteps)
            Crack_fraction = 0; // Start each dp iteration with Zero crack fraction
        } // end if max_cFaces_fraction > 0 for crack_loop

            /// Increment of special faces fraction
            dp = max_sFaces_fraction / (double) number_of_fsteps; // Face fraction increment
        } /// END of loop : for (number_of_fsteps)
    } //if max_sFaces_fraction > 0
    else { cout << "WARNING!!! max_sFaces_fraction = 0" << endl; }

        // closing all the ofstreams
        OutFLfile.close();
        OutElCondfile.close();
        OutSFile.close();

    } /// end SIMULATION MODE if (FRACTURE_LOOP)
 */