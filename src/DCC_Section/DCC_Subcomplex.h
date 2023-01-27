/// Attached user defined C++ libraries:
///-------------------------------------
#include "Planecut_Functions.h"
///-------------------------------------

//std::vector <unsigned int> DCC_Subcomplex(subcomplex new_cut, std::vector<unsigned int> const &s_faces_sequence, std::vector<unsigned int> const &c_faces_sequence) {
subcomplex DCC_Subcomplex(subcomplex &new_cut, std::vector<unsigned int> const &s_faces_sequence, std::vector<unsigned int> &sub_faces_sequence, std::vector<unsigned int> const &c_faces_sequence, double a_coeff = 0.0, double b_coeff = 0.0, double c_coeff = 1.0, double D_coeff = 0.6) {
// sub_grains_sequence - all grains in the subcomplex, sub_faces_sequence - all faces in the subcomplex, common_faces_sequence - all faces common for two grains in the subcomplex, s_sub_faces_sequence - special faces, c_sub_faces_sequence - induced (fractured, for instance) faces
std::vector <unsigned int>  sub_grains_sequence, common_faces_sequence, s_sub_faces_sequence, c_sub_faces_sequence;

    SpMat GFS = SMatrixReader(paths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains
    SpMat AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    ///  Full symmetric AGS matrix instead of triagonal
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose()));

    /// Vertex coordinates reader into triplet double vector

    vector<tuple<double, double, double>> vertex_coordinates = vertex_coordinates_vector;
    vector<tuple<double, double, double>> subcomplex_grain_coordinates, subcomplex_face_coordinates;
    vector<tuple<double, double, double>> common_faces_coordinates;

/// All subcomplex grains (subcomplex_grain_sequence)
/// Grains in a plane
    sub_grains_sequence = DCC_Plane_cut(a_coeff, b_coeff, c_coeff, D_coeff);
//REPAIR    for (auto u : sub_grains_sequence)  cout << "sub_grains_sequence_grains: " << u << endl; //    cout << "grain_coordinates_vector.size(): " << grain_coordinates_vector.size() << endl;

/// Common grain coordinates
//--------------------------------------------
if (sub_grains_sequence.size() > 0) {
    for (auto subgc: sub_grains_sequence)
        subcomplex_grain_coordinates.push_back(grain_coordinates_vector.at(subgc));
}
else cout << "Caution! sub_grains_sequence.size() = 0 in DCC_Subcomplex.h" << endl;

/// All subcomplex faces (sub_faces_sequence)
    sub_faces_sequence.clear();
    common_faces_sequence.clear();

if (sub_grains_sequence.size() > 0) {
    unsigned int face_counter = 0;
    for (unsigned int l = 0; l < CellNumbs.at(2); l++) { // for each GB
        for (auto grain_id : sub_grains_sequence) {
            if (GFS.coeff(l, grain_id) == 1) {
                sub_faces_sequence.push_back(l);
                ++face_counter;
            }
        } // end for (auto grain_id : sub_grains_sequence)
//      if (face_counter > 1) common_faces_sequence.push_back(l); face_counter = 0;
    } // end of for (unsigned int l = 0; l < CellNumbs.at(2); l++)
} // end of if(sub_grains_sequence.size() > 0)


    cout << "sub_grains_sequence size " << sub_grains_sequence.size() << endl;
    if (sub_grains_sequence.size() > 0) {
//        unsigned int face_counter = 0;
        sub_faces_sequence.clear();
        s_sub_faces_sequence.clear();
            for (auto grain_id : sub_grains_sequence) {
              ///  for (auto grain_id2 : sub_grains_sequence)
                   /// if (AGS.coeff(grain_id1, grain_id2) == 1 && grain_id1 != grain_id2) {
                     ///   cout << " grain_id1 " << grain_id1 << " grain_id2 " << grain_id2 << endl;

                        /*
                       grain3D grain1(grain_id1); grain3D grain2(grain_id2);
                        grain1.Set_GBs_list(grain_id1, GFS);
                        grain2.Set_GBs_list(grain_id2, GFS);
                        cout << " New grains: " << endl;
                        cout << " grain_id1 " << grain_id1 << " grain_id2 " << grain_id2 << endl;
                        cout << " grain_id1: "<< endl;
                        for (auto gr1 : grain1.Get_GBs_list())
                        cout << gr1 << endl;
                        cout << " grain_id2: "<< endl;
                        for (auto gr2 : grain2.Get_GBs_list())
                            cout << gr2 << endl;
*/
                        for (unsigned int l = 0; l < CellNumbs.at(2); l++) {
                            if (GFS.coeff(l, grain_id) == 1) {
                                sub_faces_sequence.push_back(l);
                                if (find(s_faces_sequence.begin(), s_faces_sequence.end(), l) == s_faces_sequence.end())
                                    s_sub_faces_sequence.push_back(l);
                            } //end of if (GFS.coeff(l, grain_id) == 1)
                            }

                    } // end of for (auto grain_id : sub_grains_sequence)
    } // end of if (sub_grains_sequence.size() > 0)

cout << "subcomplex faces sequence size: " << sub_faces_sequence.size() << endl;
cout << "s_sub_faces_sequence size: " << s_sub_faces_sequence.size() << endl;

/*
    for (auto face_id: sub_faces_sequence)
        if (count(sub_faces_sequence.begin(), sub_faces_sequence.begin() + sub_faces_sequence.size(), face_id) > 1)
            common_faces_sequence.push_back(face_id);
*/
//cout << "common faces sequence size: " << common_faces_sequence.size() << endl;


/// Common subcomplex faces (common_faces_sequence)
/*
if (sub_faces_sequence.size() > 0) {
    bool is_common = 0;
    unsigned int face_counter = 0;

    for (auto face_id: sub_faces_sequence) {
        is_common = 0;
        cout << "faces_sequence: " << count(sub_faces_sequence.begin(), sub_faces_sequence.end(), face_id) << endl;

        if (count(sub_faces_sequence.begin(), sub_faces_sequence.end(), face_id) > 1) {
            common_faces_sequence.push_back(face_id);
        }

        face_counter = face_id;
        for (auto face_id2: sub_faces_sequence) {
            if(face_id2 == face_counter) ++face_counter;
            if(face_counter > 2) is_common = 1;
        }
        if(is_common) cout << "Is common face here !" << endl;
    }
}
 */

    /// Common face coordinates
    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber)
        if(std::find(common_faces_sequence.begin(), common_faces_sequence.end(), fnumber) != common_faces_sequence.end())
            common_faces_coordinates.push_back(find_aGBseed(fnumber, paths, CellNumbs, grain_coordinates_vector)); // tuple<double, double, double> NewSeed_coordinates = make_tuple(0, 0, 0);

//    for(auto cfs : common_faces_sequence) cout << "common_faces_sequence : " << cfs << endl;
    cout << "common faces coordinates size: " << common_faces_coordinates.size() << endl;

    /// Special faces
    if (s_faces_sequence.size() > 0)
        for (unsigned int face_number : s_faces_sequence)
            if (std::find(sub_faces_sequence.begin(), sub_faces_sequence.end(), face_number) != sub_faces_sequence.end())
                s_sub_faces_sequence.push_back(face_number);
    cout << "subcomplex special faces size: " << s_sub_faces_sequence.size() << endl;


    /// Induced faces
    if (c_faces_sequence.size() > 0)
        for (unsigned int face_number : c_faces_sequence)
            if (std::find(sub_faces_sequence.begin(), sub_faces_sequence.end(), face_number) != sub_faces_sequence.end())
                c_sub_faces_sequence.push_back(face_number);

    /// Setting all quantities to the subcomplex new_subPCC with id = 0
    new_cut.Set_grains_sequence(sub_grains_sequence);
    new_cut.Set_faces_sequence(sub_faces_sequence);
    new_cut.Set_common_faces_coordinates(common_faces_coordinates);
    new_cut.Set_sub_grain_coordinates(subcomplex_grain_coordinates);
    new_cut.Set_sfaces_sequence(s_sub_faces_sequence); //special faces
    new_cut.Set_cfaces_sequence(c_sub_faces_sequence); //cracked (induced) faces

return new_cut;
}

/// Heap ///
//std::vector<unsigned int> const   &c_faces_sequence