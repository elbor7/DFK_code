///================================ A part of the DCC Subcomplex module =============================================================///
///=================================================================================================================================///
/** The library contains functions providing k-Cells IDs for a given subcomplex cut rules               **/
///================================================================================================================================///

/// The [parallelised] function output grain_ids from the plain cut (a,b,c,D)
std::vector<unsigned int> DCC_Plane_cut (double a_coeff, double b_coeff, double c_coeff, double D_coeff) {
/// The plane parameters: a_coeff*X + b_coeff*Y + c_coeff*Z = D
    vector<unsigned int> planecut_grains;

    // Obtaining Faces (coloumns) - Edges (rows) incidence matrix from file paths.at(5 + (dim - 3))
    SpMat ENS = SMatrixReader(paths.at(4 + (dim - 3)), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
    //SpMat AES = SMatrixReader(paths.at(1 + (dim - 3)), (CellNumbs.at(1)), (CellNumbs.at(1))); //all Edges
    ///  Full symmetric AES matrix instead of triagonal
    //AES = 0.5 * (AES + SparseMatrix<double>(AES.transpose()));

    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1), CellNumbs.at(2)); // Edges-Faces
    //SpMat AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    ///  Full symmetric AFS matrix instead of triagonal
    //AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose()));

    SpMat GFS = SMatrixReader(paths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains
    SpMat AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    ///  Full symmetric AGS matrix instead of triagonal
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose()));

    /// Grains
    vector<grain3D> grains_list(CellNumbs.at(3),0); // vector of grains (class grain3D)

#pragma omp parallel for // parallel execution by OpenMP
    for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
        grain3D new_grain = grain3D(m);
//        cout << new_grain.grain_id << endl;
        grains_list.at(m) = new_grain;
//        cout <<  "grains_list: " << grains_list.size() << endl;
//        cout << GFS.cols() << "  " << CellNumbs.at(3) << endl;
        grains_list.at(m).Set_node_ids(m, GFS, FES, ENS);
//        cout << grains_list.at(m).grain_id << endl;
//        cout << "m = " << m << " grain list: " << grains_list.at(m).Get_node_ids(m).size() << endl;
        grains_list.at(m).Set_node_coordinates(m);
//                cout << grains_list.at(m).grain_id << endl;
    } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) - Grains

    /// For each grain minmax_coord vector grain_coordinate_extremums of two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}
#pragma omp parallel for // parallel execution by OpenMP
    for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
        vector<tuple<double, double, double>> grain_coordinate_extremums = grains_list.at(m).Get_minmax_node_coordinates(m); // get<0>(grain_coordinate_extremums.at(0)) -> MIN X coordinate, get<1>(grain_coordinate_extremums.at(0)) -> MIN Y coordinate, get<0>(grain_coordinate_extremums.at(1)) -> MAX X coordinate, etc..
//REPAIR cout <<" gain # : " << m << " " << grains_list.at(m).grain_id << endl;
//REPAIR cout << "Xmin " << get<0>(grain_coordinate_extremums.at(0)) << " Ymin " <<get<1>(grain_coordinate_extremums.at(0)) << " Zmin " << get<2>(grain_coordinate_extremums.at(0)) << endl;
//REPAIR cout << "Xmax " << get<0>(grain_coordinate_extremums.at(1)) << " Ymax " <<get<1>(grain_coordinate_extremums.at(1)) << " Zmax " << get<2>(grain_coordinate_extremums.at(1)) << endl;
        ///MinMax condition: / a_coeff*X + b_coeff*Y + c_coeff*Z = D
        if(a_coeff*get<0>(grain_coordinate_extremums.at(0)) + b_coeff*get<1>(grain_coordinate_extremums.at(0)) + c_coeff*get<2>(grain_coordinate_extremums.at(0)) < D_coeff  &&
           a_coeff*get<0>(grain_coordinate_extremums.at(1)) + b_coeff*get<1>(grain_coordinate_extremums.at(1)) + c_coeff*get<2>(grain_coordinate_extremums.at(1)) > D_coeff ) // simultaneously: z_min < D < z_max
                planecut_grains.push_back(m);
    } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
//a,b,c,d

//REPAIR    for (auto gid : planecut_grains)         cout << gid << endl;

    return planecut_grains;
}
