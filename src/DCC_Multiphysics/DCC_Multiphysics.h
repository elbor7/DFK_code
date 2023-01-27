///================================ DCC Multiphysics module ===================================================================================///
///========================================================================================================================================///
/** The function in this library assign scalar and vector quantities to PCC k-cells associated with some physical parameters.             ///
*  In particular, elastic energies, temperature, conductivity, etc.                                                          **/         ///
///=====================================================================================================================================///

/// Attached user-defined C++ libraries:
///-------------------------------------
#include "Multiphysics_Functions.h"
///-------------------------------------

using namespace std;
//using namespace Eigen;

/// User-defined function definitions:
/* No one to date */

//std::vector <double> DCC_Multiphysics(std::vector <macrocrack> &large_cracks_vector, double external_vonMizes_stress) { // std::vector<unsigned int> &s_faces_sequence, string K_type
std::vector <double> DCC_Multiphysics(macrocrack &large_crack, double external_vonMizes_stress) { // std::vector<unsigned int> &s_faces_sequence, string K_type
    std::vector <double> new_face_energies(CellNumbs.at(2),0.0); // initial values 0.0

    int set_crack_mode = 1;
    double Puasson_coeff = 0.3;
    /// Set global sample size vector here (!)
    double average_grain_size = 500.0*pow(10,-9); // [metres] | 500 nm grains (!)
    int grains_in_a_row = round(pow(CellNumbs.at(3),0.3333)); /// (!) currently the same for 3 directions
    sample_size_vector = {grains_in_a_row*average_grain_size, grains_in_a_row*average_grain_size, grains_in_a_row*average_grain_size}; // in meters!

        ///setting crack lengths (!)
        large_crack.Set_real_crack_length(sample_size_vector.at(1)); //Set_real_crack_length(double sample_size)
       double new_crack_length = large_crack.Get_real_crack_length();
//REPAIR
cout << "grains_in_a_row: " << grains_in_a_row << " sample_size_vector(0): " << sample_size_vector[0] << " macrocrack length: " << new_crack_length << endl;
        /// Call the main function
        crack_modes_stress_field(new_face_energies, set_crack_mode, new_crack_length, external_vonMizes_stress, Puasson_coeff);
    /// Surface energy of a macrocrack
        double sface_energy_matrix = 2.0, adhesion_energy_rGO = 1.0; /// crack real energies  here (!)

//     large_cracks_vector.at(itr).surface_energy = 2.0*new_crack_length*sample_size_vector.at(0)*sface_energy_matrix; // 2 (two surfaces)*L_y*Size_x
        large_crack.surface_energy = 2.0*new_crack_length*sample_size_vector.at(0)*sface_energy_matrix; // 2 (two surfaces)*L_y*Size_x
        double bridging_coeff = 4.0; /// show how many GBs contributes to the bridging effect

        double total_in_crack_sfaces_area = 0.0;
        for (unsigned int gb : large_crack.Get_sfaces_sequence())
             total_in_crack_sfaces_area  += face_areas_vector.at(gb); // total area of all sfaces
         large_crack.bridging_energy = bridging_coeff*adhesion_energy_rGO*total_in_crack_sfaces_area*sample_size_vector.at(0)*sample_size_vector.at(1);
             //large_cracks_vector.at(itr).bridging_energy = bridging_coeff*adhesion_energy_rGO*total_in_crack_sfaces_area*sample_size_vector.at(0)*sample_size_vector.at(1);

//REPAIR    cout << " macro_crack.bridging_energy: " << large_cracks_vector.at(0).bridging_energy << " macro_crack.surface_energy: " << large_cracks_vector.at(0).surface_energy <<endl; ++itr;  } // end of for( auto macro_crack : large_cracks_vector)

    return new_face_energies;
};


/*
    /// All grain coordinates
    string GCpath_string = input_folder + "grain_seeds.txt"s;
    char* GCpath = const_cast<char*>(GCpath_string.c_str());
    vector<tuple<double, double, double>> grain_coordinates = TuplesReader(GCpath);

    /// All face coordinates
    vector<tuple<double, double, double>> face_coordinates;
    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber)
            face_coordinates.push_back(find_aGBseed(fnumber, paths, CellNumbs, grain_coordinates));
*/
