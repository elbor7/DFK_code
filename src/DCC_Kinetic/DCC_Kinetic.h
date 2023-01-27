///================================ DCC Kinetic module =============================================================///
///=======================================================================================================================///
/** The functions in this library generate different quasi-random kinetic processing the elements of the pre-constructed ///
*   polyhedral cell complex (PCC)                                                                                         **/
///=====================================================================================================================///
/// 'W' for the 3D one-layer film, 'P' for the Ising-like model of Plasticity, 'F' for the Ising-like model of Fracture
using namespace std;
using namespace Eigen;

// Triplets in the form T = T(i,j,value), where i and j element's indices in the corresponding dense matrix
typedef Triplet<double> Tr; // <Eigen library class> Declares a triplet's type name - Tr
typedef SparseMatrix<double> SpMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles name - SpMat
typedef tuple<double, double, double> Tup; // Eigen library class
typedef pair<double, double> Pr; // Eigen library class

/// Attached user defined C++ libraries:
///-------------------------------------
#include "Kinetic_Functions.h"
///-------------------------------------

std::vector <unsigned int> DCC_Kinetic(std::vector<double> &face_elastic_energies, std::vector<unsigned int> &s_faces_sequence, macrocrack &large_crack, string K_type) {

/// Output vector of the "very special" faces as the result of the module work
    std::vector <unsigned int> face_sequence;

    /// Type of the Kinetic tool from config.txt
    char* ktype = const_cast<char*>(K_type.c_str()); // type of the kinetic function from config.txt
    char* odir = const_cast<char*>(output_folder.c_str()); // const_cast for output directory
    char cktype = *ktype;

/// Declaration of FUNCTIONS, see the function bodies at the end of file.
    Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols); // The function read any matrices from lists of triplets and create corresponding sparse matrices

//    std::vector<unsigned int> VectorReader(char* FilePath); // The function read any matrices from lists of triplets and create corresponding sparse matrices
//    std::vector<int> confCout(char* config, vector<int> const& configuration); // Read and Output configuration

/// Reading vector from the file "number_of_cells" the numbers of cells od different types ///
// ::      vector components: [0] - Nodes, [1] - Edges, [2] - Faces, [3] - Grains ::     ///
    // File path with the amounts of the different cells  (1st line for Nodes, 2nd line for Edges, 3rd line for Faces and (in 3D) 4th line for grains)
/*    string  ncells = input_folder + "number_of_cells.txt"s; char* number_of_cells = const_cast<char*>(ncells.c_str());
    std::vector<unsigned int> CellNumbs = VectorReader(number_of_cells);
    //Screen output for the numbers of cells in the DCC
    cout << "The number of different k-cells in the DCC:" << endl;
  */
////////////////////// Matrices initialisation part //////////////////////
    std::vector<Tr> SFaces_Triplet_list; // Probe vector of triplets

////////==================== Reading from files of all the adjacency matrices of the DCC ===================////
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AG - Grains (Volumes or 3-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MGF - Grains - Faces (3-Cells to 2-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values // odir - source path
/// Adjacency sparse matrix for nodes /// Adjacency sparse matrix for edges /// Adjacency sparse matrix for faces /// Adjacency sparse matrix for grains
//    SpMat ANS(CellNumbs.at(0), CellNumbs.at(0)), AES(CellNumbs.at(1), CellNumbs.at(1)),
//            AFS(CellNumbs.at(2), CellNumbs.at(2)), AGS(CellNumbs.at(3), CellNumbs.at(3));
    SpMat AFS(CellNumbs.at(2), CellNumbs.at(2));
//    ANS = SMatrixReader(paths.at(_2D_), (CellNumbs.at(0)), (CellNumbs.at(0))); //all Nodes
//    ANS = 0.5 * (ANS + SparseMatrix<double>(ANS.transpose())); // Full matrix instead of triagonal
//    AES = SMatrixReader(paths.at(1 + (dim - 3)), (CellNumbs.at(1)), (CellNumbs.at(1))); //all Edges
//    AES = 0.5 * (AES + SparseMatrix<double>(AES.transpose())); // Full matrix instead of triagonal
    AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose())); // Full matrix instead of triagonal
//   AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
//    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal
/// Incidence sparse matrix for Edges and Nodes /// Incidence sparse matrix for Faces and Edges /// Incidence sparse matrix for Grains and Faces
//    SpMat ENS(CellNumbs.at(0), CellNumbs.at(1)), FES(CellNumbs.at(1), CellNumbs.at(2)),
//            GFS(CellNumbs.at(2),CellNumbs.at(3));
//    ENS = SMatrixReader(paths.at(4 + (dim - 3)), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
    SpMat   FES = SMatrixReader(paths.at(5 + (dim - 3)), (CellNumbs.at(1)), (CellNumbs.at(2))); //all Edges-Faces
//    GFS = SMatrixReader(paths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains
 /// =================================================================================================///
////////////////////////////////////// KINETIC PROCESSES ///////////////////////////////////////////////

    /// #2# Kinetic type FRACTURE (CRACKING)
    if (cktype == 'F' ) { // Fracture
        //face_sequence = DCC_Kinematic_cracking(s_faces_sequence, AFS, FES);
        face_sequence = DCC_Kinetic_cracking(face_elastic_energies, s_faces_sequence, large_crack, AFS, FES);
    } ///End of 'Fracture' type simulations

    /// #2# Kinetic type WEAR
    else if (cktype == 'W') { // Wear
        double  ShearStress = 3.0*pow(10,8);
        /// Creation/Reading of grain orientations
        double Ori_angle = 0;
        vector<Tup> Grain_Orientations; // Vectors of triplets (in 2D {x,y,0}) for grain orientations
        for (unsigned int k = 0; k < CellNumbs.at(3); k++) {
            Ori_angle = rand() % 180; // Random generation of angle
            //if (Ori_angle < 90)
            Grain_Orientations.push_back(make_tuple(cos(Ori_angle),sin(Ori_angle),0));
            // cout << get<0>(Grain_Orientations.at(k)) << "\t" << get<1>(Grain_Orientations.at(k)) << "\t" <<get<2>(Grain_Orientations.at(k)) << "\t" << endl;
        }
        /// ======= Kinetic function ===================>>
//        DCC_Kinetic_Wear(ShearStress, Grain_Orientations, FES, CellNumbs, input_folder, odir);
        DCC_Kinetic_Wear( ShearStress, Grain_Orientations, FES);

        /// ===== Grain orientations output =========>
        ofstream GrainOrientationsStream;

        GrainOrientationsStream.open(odir + "GrainOrientations.txt"s, ios::trunc);
        if (GrainOrientationsStream) {
            GrainOrientationsStream << "Grain Orientations" << endl;
            for (auto go : Grain_Orientations)
                GrainOrientationsStream << get<0>(go) << "\t" << get<1>(go) << endl; //<< "\t" <<get<2>(go)
            GrainOrientationsStream.close();
        } else cout << "Error: No such a directory for\t" << odir + "GrainOrientations.txt"s << endl;

    } ///End of 'Wear' type simulations
/// -----------------------------------------------------------------------------------------------------------------///

    /// #3# Kinetic type PLASTICITY
    else if (cktype == 'P') { // Plasticity
        vector<Tup> fraction_stress_temperature;

        /// ========= DCC_Kinetic main function call ==========>>
//        fraction_stress_temperature = DCC_Kinetic_Plasticity(FES, CellNumbs, input_folder, odir);
        fraction_stress_temperature = DCC_Kinetic_Plasticity( FES );

        /// Analysis and output to file
        ofstream FSTStream;
        FSTStream.open(odir + "fraction_stress_temperature.txt"s, ios::trunc);
        if (FSTStream) {
            FSTStream << "(1) Plastic strain [%] (2) Yield strength [MPa] (3) Temperature [K])" << endl;
            for (auto go : fraction_stress_temperature)
                FSTStream << get<0>(go) * 100.0 << " \t" << get<1>(go)/pow(10,6) << " \t" << get<2>(go) << endl;
            FSTStream.close();

        } else cout << "Error: No such a directory for\t" << odir + "fraction_stress_temperature.txt"s << endl;
//        for (auto lk : fraction_stress_temperature) if (get<0>(lk)*CellNumbs.at(2)*Burgv/DCC_size > 0.002) { cout << "Nano-slips fraction\t" << get<0>(lk) << "\tYield strength\t" << get<1>(lk) << "\tTemperature\t" << get<2>(lk) << endl;              break; }

    } ///End of 'Plasticity' type simulations
/// -----------------------------------------------------------------------------------------------------------------///
    else cout << "ERROR [DCC_Kinetic] : unknown simulation type - please replace with 'W or P' " << endl; return face_sequence;

// Closing and deleting

    return face_sequence;
} /// The end of DCC_Kinetic()

/// ================================== Related functions ==================================///
