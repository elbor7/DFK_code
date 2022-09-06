///******************************************************************************************************************************///
///*                                                                                                                            *///
///******************************************************************************************************************************///
///********************************** Dr Elijah Borodin, Manchester, UK, Spring 2022   ***********************************************///
///******************************************************************************************************************************///
///****************************************   DCCAnalyser(c) utility   ********************************************************///
///******************************************************************************************************************************///

/// Standard (STL) C++ libraries:
///------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <tuple>
#include<numeric>
#include <algorithm>
///-------------------------------
/// Attached user defined C++ libraries:
///------------------------------------------
// #include <Eigen/Core>
// #include <Eigen/Dense>
#include <Eigen/SparseCore> // Eigen source: https://eigen.tuxfamily.org/ (2022)
#include <Spectra/MatOp/SparseGenMatProd.h> // Spectra source: https://spectralib.org/ (2022)
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
///------------------------------------------
using namespace std; //Standard namespace
using namespace Eigen;
using namespace Spectra;

// Triplets in the form T = T(i,j,value), where i and j element's indices in the corresponding dense matrix
typedef Triplet<double> Tr; // <Eigen library class> Declares a triplet's type name - Tr
typedef SparseMatrix<double> SpMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles name - SpMat

/// Declaration of Global variables
int dim, number_of_types;
std::vector<unsigned int> CellNumbs;
vector<char*> paths;
double max_sFaces_fraction, max_cFaces_fraction;
string input_folder, output_folder; // input and output folders from file config.txt file
ofstream OutFLfile, OutElCondfile, OutSFile; // ofstreams for data output


//vector<bool> SChar_config; // Characterisation module configuration

/* Various useful functions (must be here first in the list)*/
#include "src/DCC_SupportFunctions.h"
/* Processing module assigned new id's for the various DCC elements (Faces, Edges,...) */
/* As its output module generates sequences of element's id's and the State vector of Faces                  */
#include "src/DCCProcessing/DCCProcessing.h"
/* Kinetic module assign some values for the scalar or vector variables defined on the DCC elements (Faces, Edges,...) */
/* As its output module generates one or several sequences of element's id's and the State vector of Faces                  */
#include "src/DCCKinetic/DCCKinetic.h"
/* Characterisation module counts various structural characteristics of special substructures defined on the DCC elements */
#include "src/DCCCharacterisation/StructureCharacterisation.h"

/// Declaration of FUNCTIONS, please see the function bodies at the end of file.
std::vector<double> confCount(char* config, string &Processing_type, string &Kinetic_type, string &input_folder, string &output_folder); // Read and output the initial configuration from the config.txt file
string SIMULATION_MODE(char* config); // Check the SIMULATION MODE type in the config.txt file
bool ProcessingON(char* config, bool time_step_one); // Check the Processing module status (On/OFF) in the config.txt file
bool CharacterisationON(char* config, bool time_step_one); // Check the Structure Characterisation module status (On/OFF) in the config.txt file
bool KineticON(char* config, bool time_step_one); // Check the Kinetic module status (On/OFF) in the config.txt file
void eraseSubStr(std::string & mainStr, const std::string & toErase); //support function: erase the first occurrence of given substring from main string

///*........................................................................................    Main    ................................................................*///
int main() {
    cout << "-------------------------------------------------------------------------" << endl;

/// File path to the configuration profile
// Still not working with Windows:(
    using std::__fs::filesystem::current_path; //to obtain current working directory
    string MainPath = current_path();
    eraseSubStr(MainPath, "cmake-build-debug"s);
    eraseSubStr(MainPath, "buildtree"s);
    string  config = MainPath + "config.txt"s; char* confpath = const_cast<char*>(config.c_str());

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
    // The source directory and simulation type from file config.txt
    string P_type; // 'R', 'S', 'D', 'I' or 'E' :: This char define the process type: 'R' for Random, 'S' for maximum configuration Entropy production, 'D' for DDRX process (DCC retessellations with the new seeds), 'I' for the Index-based generation modes, 'E' for experimental data obtained by EBSD
    string K_type; // 'W', 'P' or 'F' :: This char define the Kinetic process type: // 'W' for the 3D one-layer film, 'P' for the Ising-like model of Plasticity, 'F' for the Ising-like model of Fracture
    vector<double> ConfigVector = confCount(confpath, P_type, K_type, input_folder, output_folder);
    char* indir = const_cast<char*>(input_folder.c_str()); // const_cast for input directory
//    char* odir = const_cast<char*>(output_folder.c_str()); // const_cast for output directory
    dim = ConfigVector.at(0); // Space dimension of the problem (dim = 2 or 3);
/// Below the file names with the sparse DCC matrices which must already exit in the input_folder and have the same names!
    string ssd0 = input_folder + "A0.txt"s, ssd1 = input_folder + "A1.txt"s, ssd2 = input_folder + "A2.txt"s, ssd3 = input_folder + "A3.txt"s, ssd4 = input_folder + "B1.txt"s, ssd5 = input_folder + "B2.txt"s, ssd6 = input_folder + "B3.txt"s,
            seeds = input_folder + "seeds.txt"s, NewSeeds = input_folder + "NewSeeds/NewSeeds.txt"s;
//The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
    paths.push_back(const_cast<char*>(ssd0.c_str())); paths.push_back(const_cast<char*>(ssd1.c_str())); paths.push_back(const_cast<char*>(ssd2.c_str())); paths.push_back(const_cast<char*>(ssd3.c_str()));
    paths.push_back(const_cast<char*>(ssd4.c_str())); paths.push_back(const_cast<char*>(ssd5.c_str())); paths.push_back(const_cast<char*>(ssd6.c_str()));
    paths.push_back(const_cast<char*>(seeds.c_str())); paths.push_back(const_cast<char*>(NewSeeds.c_str()));
// File path with the amounts of the different cells  (1st line for Nodes, 2nd line for Edges, 3rd line for Faces and (in 3D) 4th line for grains)
    string  ncells = input_folder + "number_of_cells.txt"s; char* number_of_cells = const_cast<char*>(ncells.c_str());
    bool time_step_one = 1; // id for the first step of iteration MAX fraction of Faces
/// Output paths.vector to console out
    int npath = 0;
    cout << "___________________________________________________________________________" << endl;
    for (auto m : paths) cout <<"[" << npath++ << "]" << " paths:\t" << m << endl;

/// Principal variables
///_____________________________________________________________________________________
    // Read all the number of cells from file number_of_cells
    // :: CellNumbs vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Grains number
    CellNumbs = VectorReader(number_of_cells);

    /// Special feature for 2D case:
    if (dim == 2) { CellNumbs.push_back(CellNumbs.at(2)); CellNumbs.at(2) = CellNumbs.at(1); CellNumbs.at(1) = CellNumbs.at(0); CellNumbs.at(0) = NULL; }

    /// CellNumbs output::
    cout << "=========================================================================" << endl;
    unsigned int t_length = 0;
    for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl;

    // State vector | Special faces IDs
    vector <unsigned int> State_sVector(CellNumbs.at(2), 0), current_State_sVector(CellNumbs.at(2), 0), State_cVector(CellNumbs.at(2), 0), current_State_cVector(CellNumbs.at(2), 0);
    // Order of newly generated special faces | Process time
    std::vector <unsigned int> special_faces_sequence, current_sfaces_sequence, crack_faces_sequence, current_cracks_sequence; // Variable sequences (in order of their generation) of special Faces and Cracks
    // The number of special Face (2-cells) types
    number_of_types = ConfigVector.at(1); // number of different special cell id's
    // MAX fraction of special Faces | Calculation limit
    max_sFaces_fraction = ConfigVector.at(2), max_cFaces_fraction = ConfigVector.at(3);

cout << "=========================================================================" << endl << "\t\t\t\t\t[\tStart of the DCC Processing Tool\t]\t\t\t\t\t" << endl << "-------------------------------------------------------------------------" << endl;
    /// Output (ofstream) paths for different files from Characterisation module
    string output_TJs_dir = output_folder + "TJsLab_TJsTypes.txt"s; char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str());  // From string to char for the passing folder path to a function
    string output_Entropy_dir = output_folder + "TJsLab_TJsEntropy.txt"s; char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str());
    string Face_Cond_dir = output_folder + "Face_Conductivity.txt"s; char* FCond_dir = const_cast<char*>(Face_Cond_dir.c_str());
    string Face_Laplacian_dir = output_folder + "SpecialFacesLaplacian.txt"s; char* FLap_dir = const_cast<char*>(Face_Laplacian_dir.c_str());

    /// Creation of the two files with TJs fraction and Configuration entropies output
    // ofstream for entropies OutSFile
    ofstream OutSFile; OutSFile.open(cEntropy_dir, ios::trunc);
    OutSFile << "\t(1)Fraction_of_special_Faces\t" << " " << "\t(2)Fraction_of_cracks\t" << " " << "\t(3)Von-Neumann_Special_Faces_entropy\t" << " " << "\t(4)j0_TJs\t"<< " " << "\t(5)j1_TJs\t"<< " " << "\t(6)j2_TJs\t"<< " " << "\t(7)j3_TJs\t" << " " << "\t(8)p*=j1+2*j2+3*j3/3\t" << " " << "\t(9)Special_Faces_Entropy\t" << " " << "\t(10)Special_Faces_median_entropy\t" << " " << "\t(11)Special_Faces_skrew_entropy\t" << " " << "\t(12)Special_Faces_informativeness\t" << " "
             << "\t(13)Von-Neumann_Cracked_Faces_entropy\t" << " " << "\t(14)jc0_cTJs\t"<< " " << "\t(15)jc1_cTJs\t"<< " " << "\t(16)jc2_cTJs\t"<< " " << "\t(17)jc3_cTJs\t" << " " << "\t(18)f*=jc1+2*jc2+3*jc3/3\t"<< " " << "\t(19)Special_Cracks_Entropy\t" << " " << "\t(20)Cracked_Face_median_entropy\t" << " " << "\t(21)Cracked_Face_skrew_entropy\t" << " " << "\t(22)Cracked_Face_informativeness\t" << endl;
    OutSFile << endl;
    OutSFile.close();
    // ofstream for TJs OutTJsFile
    ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::trunc);
    OutTJsFile << "\t(1)Fraction_of_special_Faces\t" << " " << "\t(2)j0_TJs\t"<< " " << "\t(3)j1_TJs\t"<< " " << "\t(4)j2_TJs\t"<< " " << "\t(5)j3_TJs\t" << " " << "\t(6)Configurational_Face_Entropy\t" << " " << "\t(7)Face_Entropy_Median\t" << " " << "\t(8)Face_Entropy_Skrew\t" << endl;
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

/// ==========================================================================================================================================
/// ================================================= THE LIST MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
    if ( SIMULATION_MODE(confpath) == "SIMULATION MODE(LIST)"s) { /// SIMULATION MODE: #LIST
/// I: DCC_Processing module
        if ( ProcessingON(confpath, time_step_one)) { // if DCC_Processing is SWITCH ON in the config.txt file
//            DCC_Processing3D(State_Vector, special_faces_sequence, Processing_type, max_sFaces_fraction, number_of_types, CellNumbs, paths, odir, ip_index);
            DCC_Processing3D( State_sVector, special_faces_sequence, P_type);
        } // end if(ProcessingON(confpath))
        else {
            /// ====== Reading from file instead of DCC_Processing3D ( ) if it is off ============>
            string SFS_path = output_folder + "Smax/SpecialGrainBoundaries.txt"s;
            char* SFS_dir = const_cast<char *>(SFS_path.c_str());
            special_faces_sequence = Smax_sequence_reader(SFS_dir); // a function from Processing_Functions.h
            for (auto itr: special_faces_sequence) State_sVector.at(itr) = 1;
        } // end of else

/// III: DCC_Kinetic module
        if (KineticON(confpath, time_step_one)) { // if DCC_Kinetic is SWITCH ON in the config.txt file
//            DCC_Kinetic(Kinetic_type, special_faces_sequence, paths, indir, odir);
            DCC_Kinetic(special_faces_sequence, K_type);
        }// end if(KineticON(confpath))

/// II: DCC_Characterisation module
        if (CharacterisationON(confpath, time_step_one)) {
            cout << "START of DCC Structure Characterisation Module" << endl;
//            DCC_StructureCharacterisation(State_sVector, special_faces_sequence, crack_faces_sequence, ConfigVector, CellNumbs, paths, odir, OutFLfile, OutElCondfile, OutSFile);
            DCC_StructureCharacterisation(State_sVector, special_faces_sequence, crack_faces_sequence, ConfigVector);

        }// end if(CharacterisationON(confpath))
    } /// end SIMULATION MODE else (LIST)

/// ===========================================================================================================================================
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
    long number_of_fsteps = 300, number_of_csteps = 30; /// NUMBER OF STEPS for loops over the fractions of special Faces and Cracks
    double Face_fraction = 0.01, Crack_fraction = 0.00; /// INITIAL fractions of special Faces and Cracks
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
                crack_faces_sequence = DCC_Kinetic(current_sfaces_sequence, K_type);
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
//               DCC_StructureCharacterisation(current_State_sVector, current_sfaces_sequence, current_cracks_sequence, ConfigVector, CellNumbs, paths, odir, OutFLfile, OutElCondfile, OutSFile);
                 DCC_StructureCharacterisation(current_State_sVector, current_sfaces_sequence, current_cracks_sequence, ConfigVector);

                } else break;

                // the following time steps
                time_step_one = 0;
                /// Increment of crack fraction
                df = max_cFaces_fraction / (double) number_of_csteps;
            } /// END of loop : for (number_of_csteps)
            Crack_fraction = 0; // Start each dp iteration with Zero crack fraction
        } // end if max_cFaces_fraction > 0 for crack_loop
            else {
            /// ====== Structure Characterisation module ========>
            if (CharacterisationON(confpath, time_step_one)) {
//               DCC_StructureCharacterisation(current_State_sVector, current_sfaces_sequence, current_cracks_sequence, ConfigVector, CellNumbs, paths, odir, OutFLfile, OutElCondfile, OutSFile);
                DCC_StructureCharacterisation(current_State_sVector, current_sfaces_sequence, current_cracks_sequence,
                                              ConfigVector);

            } else break;
        }
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

/// ===== Elapsing time ================>
    unsigned int end_time = clock();
    double fulltime = (double) end_time;
    cout << dim << "D " << "Runtime is equal to  " << fulltime <<  "  seconds" << endl;

    cout << "-------------------------------------------------------------" << endl << "The end of the VoroCAnalyser program" << endl << "=============================================" << endl ;

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
                if (time_step_one == 1) cout << "DCC_Processing SWITCHED ON"s << endl;
                return isProcessingON;
            }
    } else cout << "ProcessingON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "DCC_Processing SWITCHED OFF"s << endl;
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
                if (time_step_one == 1) cout << "DCC_Kinetic SWITCHED ON"s << endl;
                return isKineticON;
            }
    } else cout << "KineticON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "DCC_Kinetic SWITCHED OFF"s << endl;
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
            if (time_step_one == 1)    cout << "DCC_StructureCharacterisation SWITCHED ON"s << endl;
                return isCharacterisationON;
            }
    } else cout << "isCharacterisationON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "DCC_StructureCharacterisation SWITCHED OFF"s << endl;
    return isCharacterisationON;
}

/// Reading and Output of the configuration file
std::vector<double> confCount(char* config, string &Processing_type, string &Kinetic_type, string &input_folder, string &output_folder) {
    vector<double> res;

    string line;
    double p_max = 0, f_max = 0;

    ifstream inConf(config);
    if (inConf) { // If the file was successfully open, then
        while (getline(inConf, line, '\n')) {
            // Number of types
            for (auto it: line) {
                if (it == '&') Processing_type = line.at(0); // simulation type
                else if (it == '`') Kinetic_type = line.at(0); // simulation Kinetic type

                else if (it == '@') res.push_back(line.at(0) - '0'); // dimension of the problem; res[1]
                else if (it == '!') res.push_back(line.at(0) - '0'); // number of special Face types; res[0]

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
    cout << "The problem dimension that is the maximum dimension k_max of k-Cells:\t | "s << res.at(0) << endl;
    cout << "Calculation type ('R', 'W', 'S', 'F', 'I' or 'E'):\t\t\t\t\t\t | "s << Processing_type << endl;
    cout << "Kinetic type ('W' or 'I'):\t\t\t\t\t\t\t\t\t\t\t\t | "s << Kinetic_type << endl;
    cout << "The number of special Face (2-cells) types:\t\t\t\t\t\t\t\t | " << res.at(1) << endl;
    cout << "MAX fraction of Faces (calculation limit): \t\t\t\t\t\t\t\t | " << res.at(2) << endl;
    cout << "MAX fraction of Cracks (calculation limit):\t\t\t\t\t\t\t\t | " << res.at(3) << endl;
    cout << endl;
    cout << "Input folder:\t" << input_folder << endl;
    cout << "Output folder:\t" << output_folder << endl;
    cout << endl;
    cout << "Nodes types statistics, indices and configuration entropy:     "; if (res.at(4) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Edges types statistics, indices and configuration entropy:     "; if (res.at(5) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Faces types statistics and structural indices:                 "; if (res.at(6) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain types statistics, indices and configuration entropy:     "; if (res.at(7) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Node Laplacian with its spectrum for the Nodes network:       "; if (res.at(8) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Edge Laplacian with its spectrum for the Edges network:       "; if (res.at(9) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Face  Laplacian with its spectrum for the Faces network:       "; if (res.at(10) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain Laplacian with its spectrum for the Grains network:      "; if (res.at(11) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Tutte polynomial for the special network:                      "; if (res.at(12) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;

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
/* Two functions (dependent on the problem's spatial dimension - 2D or 3D) are launching here with the arguments of the paths for
 * adjacency AN, AE, AF, (AG) and incidence (boundary operators) BEN, BFE, (BGF) matrices. */

/*void Manual_input () {
 do { ///Manual user input of the space dimension value (dim)
    cout << " Please, input the dimension of the problem: 3 for 3D and 2 for 2D "s << endl;
    cin >> dim;
    cout << "The problem dimension is\t" << dim << endl;
    if (dim != 2 && dim != 3) cout << "Input Error: Please retype 2 or 3" << endl;
    }while (dim != 2 && dim != 3); */
/* do { ///Manual user input of the simulation type
    cout << " Please, input the symbol of particular simulation type of the problem: E (experiment), R (rndom), S (entropy maximisation) or I (Ising model):"s << endl;
    cin >> simulation_type;
    cout <<"The simulation type is\t"<< simulation_type << endl;
    if (simulation_type != 'E' && simulation_type != 'R' && simulation_type != 'S' && simulation_type != 'I') cout << "Input Error: Please retype 'E', 'R', 'S' or 'I' for the specific simulation type" << endl;
    }while (simulation_type != 'E' && simulation_type != 'R' && simulation_type != 'S' && simulation_type != 'I');  */
/* /// Manual user input of the DCC files folder path
    cout << " Please, input the name of the folder where the DCC source files are (like 1k3cells/ ):"s << endl;
    cin >> problem_folder_path; */
/* /// Manual user input of the simulation results output folder path
    cout << " Please, input the simulation results output folder path (like test/ ):"s << endl;
    cin >> results_output_folder_path;
return;
} */