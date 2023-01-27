/// Classes of special objects related to defect structures on a PCC elements

#ifndef agglomeration_H
#define agglomeration_H
#ifndef subcomplex_H
#define subcomplex_H

/// #0# The class of Grain Boundaries in a PCC
class grain_boundary{

public:
    unsigned int GB_id; // grain boundary ID

    grain_boundary(unsigned int GB_new_id) { // class simple constructor
        GB_id = GB_new_id;
    }

/// GB state
    bool is_inclusion; // 0 - no, 1 - yes
    bool is_fractured; // 0 - no, 1 - yes
    bool is_agglomeration; // 0 - no, 1 - yes

/// Set values methods

    void Set_surface_energy(vector<double> GB_SE_vector ){
        surface_energy = GB_SE_vector.at(GB_id);
    }

    void Set_external_elastic_energy(vector<double> GB_EEE_vector) { /// EEE: external elastic energy vector
        external_elastic_energy = GB_EEE_vector.at(GB_id);
    }

    void Set_crack_interaction_energy(vector<double> GB_CIE_vector) {
        crack_interaction_energy = GB_CIE_vector.at(GB_id);
    }

    void Set_Bl_energy(vector<double> GB_BLE_vector) {
        Bl_energy = GB_BLE_vector.at(GB_id);
    }

    void Set_Cl_energy(vector<double> GB_CLE_vector){
        Cl_energy = GB_CLE_vector.at(GB_id);
    }

/// Get values methods
    double Get_surface_energy(unsigned int GB_id){
        if(surface_energy != 0)
            return surface_energy;
        else {
            Set_surface_energy(GB_SE_vector);
            return surface_energy; }
    }

    double Get_external_elastic_energy(unsigned int GB_id) { /// EEE: external elastic energy vector
        if(external_elastic_energy != 0)
            return external_elastic_energy;
        else {
            Set_external_elastic_energy(GB_EEE_vector);
            return external_elastic_energy; }
    }

    double Get_crack_interaction_energy(unsigned int GB_id) {
        if(crack_interaction_energy != 0)
            return crack_interaction_energy;
        else {
            Set_crack_interaction_energy(GB_CIE_vector);
            return crack_interaction_energy; }
    }

    double Get_Bl_energy(int GB_id){
        if(Bl_energy != 0)
            return Bl_energy;
        else {
            Set_Bl_energy(GB_BLE_vector);
            return Bl_energy; }
    }

    double Get_Cl_energy(unsigned int GB_id){
        if(Cl_energy != 0)
            return Cl_energy;
        else {
            Set_Cl_energy(GB_CLE_vector);
            return Cl_energy; }
    }

    double Get_total_energy(){

        /// total_energy; //= surface_energy + external_elastic_energy + Bl_energy + Cl_energy;
    }

private:

/// Combinatorial
    int GB_edges_number;  //equal to the number of neighbours
    int GB_nodes_number;

///Geometry
    double GB_area;
    double GB_perimeter;

    tuple<double,double,double> GB_barycentre_coordinates;

///Energies
    double surface_energy;
    double external_elastic_energy;
    double crack_interaction_energy;
    double Bl_energy;
    double Cl_energy;
    double total_energy; //= surface_energy + external_elastic_energy + Bl_energy + Cl_energy;
};



/// #2# The class of Grains in a PCC

class grain3D {
    vector<tuple<double, double, double>> minmax_node_coordinates; // a vecor containing two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}

private:
    // list of nodes
    vector<unsigned int> node_ids;

    // list of GBs
    vector<unsigned int> GBs_list;

    // list of triplets of nodes coordinated
    vector<tuple<double, double, double>> node_coordinates;

public:

    unsigned int grain_id;

    grain3D(unsigned int grain_new_id) { // constructor 1 simple
        grain_id = grain_new_id;
    }

    void Set_node_ids(unsigned int grain_id, SpMat const &GFS, SpMat const &FES, SpMat const &ENS) { // set the node ids
        /// GFS -> FES -> ENS
        if(node_ids.size() == 0) {
            for(unsigned int l = 0; l < CellNumbs.at(2); l++) {// over all Faces (l)
                if (GFS.coeff(l, grain_id) == 1) {
                    for (unsigned int j = 0; j < CellNumbs.at(1); j++) // over all Edges (j)
                        if (FES.coeff(j, l) == 1) { // at the chosen Face with ID = 'l'
                            for (unsigned int i = 0; i < CellNumbs.at(0); i++) // over all Nodes
                                if (ENS.coeff(i, j) == 1) node_ids.push_back(i); // at the chosen Face with ID = 'l'
                        } // end of if (FES.coeff(l, j) == 1)
                } // end of (GFS.coeff(m, l) == 1)
            } // end of for(unsigned int l = 0; l < CellNumbs.at(2); l++) - Faces

        }/// end of if(node_ids.size() == 0)

    }

    void Set_GBs_list(unsigned int grain_id, SpMat const &GFS) {
        for (unsigned int l = 0; l < CellNumbs.at(2); l++) // for each GB
                if (GFS.coeff(l, grain_id) == 1)
                    GBs_list.push_back(l);
    }

    vector<unsigned int> Get_GBs_list() {
            if (GBs_list.size() > 0) return GBs_list;
                else cout << "coution GBs_list.size() = 0! Please Set_GBs_list(unsigned int grain_id, SpMat const &GFS)  first!"s << endl;
            }

    /// return - vector of all node (vertices) coordinates of a grain
    void Set_node_coordinates(unsigned int grain_id) { // set the node ids from Tr = triplet list
//        for (auto  itr = node_ids.begin(); itr != node_ids.end(); ++itr)
            //if(find(node_ids.begin(), node_ids.end(), distance(node_ids.begin(), itr)) != node_ids.end())
//                node_coordinates.push_back(vertex_coordinates_vector.at(*itr)); // vector<unsigned int> node_ids;
if(node_ids.size()>0) {
//REPAIR
    for (auto ids: node_ids) {
//REPAIR        cout << "vertex_coordinates_vector size " << vertex_coordinates_vector.size() << endl;
//REPAIR        cout << "ids: " << ids << " node_coordinates size: " << get<0>(vertex_coordinates_vector.at(ids)) << endl;

        node_coordinates.push_back(vertex_coordinates_vector.at(ids)); // vector<unsigned int> node_ids;
//REPAIR cout << "grain id: " << grain_id << " node_coordinates size: " << node_coordinates.size() << endl;
    }

}
        else {
            cout << "Caution! The size of node_ids vector of the grain " << grain_id << " is 0 !" << endl;
            node_coordinates.push_back(make_tuple(0,0,0)); /// change after !!
        }

    } // the end of Set_node_coordinates

    vector<unsigned int> Get_node_ids(unsigned int grain_id) { // set the node ids
        return node_ids;
    }

    vector<tuple<double, double, double>> Get_node_coordinates(unsigned int grain_id) { // set the node ids
        if (node_coordinates.size() != 0) {
            return node_coordinates;
        } else if (node_coordinates.size() > 0) {
            Set_node_coordinates(grain_id);
            return node_coordinates;
        } else {
            throw std::invalid_argument(
                    "Please call Set_node_coordinates(unsigned int grain_id, vector<tuple<double, double, double>> const &vertex_coordinates) method first!");
            return node_coordinates;
        }
    }// end of  Get_node_coordinates() method

    /// return - vector with two tuples : { x_min, y_min, z_min; x_max, y_max, z_max} of a grain witn number grain_id
    vector<tuple<double, double, double>> Get_minmax_node_coordinates(unsigned int grain_id) { // min and max {x,y,z} values of vertices for a grain
        vector<tuple<double, double, double>> minmax_tuple;
        ///Get_node_coordinates(grain_id)
        vector<tuple<double, double, double>> tup_node_coordinates = Get_node_coordinates(grain_id); // class Grain3D function Get_node_coordinates(grain_id)
//REPAIR        cout << "Xtup_node_coordinates " << get<0>(tup_node_coordinates.at(0)) << " Ytup_node_coordinates " <<get<1>(tup_node_coordinates.at(0)) << " Ztup_node_coordinates " << get<2>(tup_node_coordinates.at(0)) << endl;

        // separating in three parts
        vector<double> x_node_coordinates, y_node_coordinates, z_node_coordinates;
        for (auto itr = tup_node_coordinates.begin(); itr != tup_node_coordinates.end(); ++itr) {
            x_node_coordinates.push_back(get<0>(*itr));
            y_node_coordinates.push_back(get<1>(*itr));
            z_node_coordinates.push_back(get<2>(*itr));
        }
//REPAIR        cout << "x_node_coordinates " << (x_node_coordinates.at(0)) << " y_node_coordinates " << (y_node_coordinates.at(0)) << " z_node_coordinates " << (z_node_coordinates.at(0)) << endl;
//        std::for_each(tup_node_coordinates.begin(), tup_node_coordinates.end(), );

// Return iterators (!)
        auto itx = std::minmax_element(x_node_coordinates.begin(), x_node_coordinates.end());
        double xmin = *itx.first;
        double xmax = *itx.second;

        auto ity = std::minmax_element(y_node_coordinates.begin(), y_node_coordinates.end());
        double ymin = *ity.first;
        double ymax = *ity.second;

        auto itz = std::minmax_element(z_node_coordinates.begin(), z_node_coordinates.end());
        double zmin = *itz.first;
        double zmax = *itz.second;

        minmax_tuple = {make_tuple(xmin, ymin, zmin), make_tuple(xmax, ymax, zmax)};
//REPAIR !!        cout << "Xmin " << get<0>(minmax_tuple.at(0)) << " Ymin " <<get<1>(minmax_tuple.at(0)) << " Zmin " << get<2>(minmax_tuple.at(0)) << endl;

        return minmax_tuple;
    }

}; // end of class grain3D

/// #3# The class of Agglomeration of defects in one special face element of a PCC
class agglomeration {
    double adhesion_energy = 0;
    double surface_energy = 0;

private:
    string atype; // like "rgo"
    unsigned int aface_number = 0;
    int          apower = 0;
    int          a_average_strip_length = 0;

public:
    agglomeration(unsigned int AFace) { // constructor 1 simple
        aface_number = AFace;
    }

    agglomeration(unsigned int AFace, int AglPower) { // constructor 2 complex
        aface_number = AFace;
        apower = AglPower;
    }

    void SetAFace(unsigned int AFace)
    {
        aface_number = AFace;
    }

    void SetAType(std::string type)
    {
        atype = type;
        if (atype == "rgo") {
            surface_energy = 0.2; // units [J/m^2]
            adhesion_energy = 0.4; // units [J/m^2]
        }

    }

    void SetAPower(vector<vector<int>> const &RW_series_vector) /// Agglomerations power
    {
        int AglPower = 0;
        //vector<int> aggl_vector(CellNumbs.at(2), 0); // state vector for agglomerations (with the size # cells initially filled with all 0s) : contains # of faces with agglomerations
        for (auto RWsfv : RW_series_vector)
            for (auto RWsf: RWsfv)
                if(RWsf == aface_number) AglPower += 1;

        apower = AglPower;
    }

    void SetAvLength(vector<vector<int>> const &RW_series_vector) /// Average length of strips related to this agglomeration
    {
        int ATotalLength = 0;
        for (auto RWsfv : RW_series_vector)
            for (auto RWsf: RWsfv)
                if(RWsf == aface_number) ATotalLength += RWsfv.size();

        a_average_strip_length = ATotalLength/ (double) GetAPower(RW_series_vector);
    }

    unsigned int GetAFace()
    {
        return aface_number;
    }

    int GetAPower()
    {
        return apower;
    }

    int GetAPower(vector<vector<int>> const &RW_series_vector)
    {
        if (apower != 0) {
            return apower;
        } else {
            int AglPower = 0;
            for (auto RWsfv : RW_series_vector)
                for (auto RWsf: RWsfv)
                    if(RWsf == aface_number) AglPower += 1;

            apower = AglPower;
            return apower;
        }
    }

            int GetAvLength()
    {
        return a_average_strip_length;
    }

    int GetAvLength(vector<vector<int>> const &RW_series_vector) {
        if (a_average_strip_length != 0) {
            return a_average_strip_length;
        } else {
            int ATotalLength = 0;
            for (auto RWsfv: RW_series_vector)
                for (auto RWsf: RWsfv)
                    if (RWsf == aface_number) ATotalLength += RWsfv.size();

                a_average_strip_length = ATotalLength / (double) GetAPower(RW_series_vector);
            return a_average_strip_length;
        }
    }

}; // end of class agglomeration

/// #4# The class of stress concentrators related to defects in one special face element of a PCC

class face_concentrator {
    double total_elastic_energy = 0;

private:
    //string conc_type; // bl or cl
    unsigned int conc_face_number = 0;
    unsigned int bl_index = 0;
    unsigned int cl_index = 0;
    double Bl_elastic_energy = 0;
    double Cl_elastic_energy = 0;
public:

    face_concentrator (unsigned int ConcFace) { // constructor 1 simple
        conc_face_number = ConcFace;
    }

    face_concentrator (unsigned int ConcFace, unsigned int Bl, unsigned int Cl) { // constructor 2 complex
        conc_face_number = ConcFace;
        bl_index = Bl;
        cl_index = Cl;
    }

    unsigned int GetBl_index()
    {
        return bl_index;
    }

    int GetCl_index()
    {
        return cl_index;
    }
}; // end of class



/// #5# The class of a SUBCOMPLEX
class subcomplex {

private:
    /// 1. Combinatorics
    std::vector <unsigned int> sub_grains_sequence;
    std::vector <unsigned int> sub_faces_sequence;
    std::vector <unsigned int> sub_nodes_sequence;
    std::vector <unsigned int> common_faces_sequence;
    std::vector <unsigned int> s_sub_faces_sequence;
    std::vector <unsigned int> c_sub_faces_sequence;

    /// 2. Geometry

    //vector<tuple<double, double, double>> vertex_coordinates;
    vector<tuple<double, double, double>> common_faces_coordinates;
    vector<tuple<double, double, double>> sub_grain_coordinates;

public:
    unsigned int subcomplex_id;
    double sub_length;
    double a_n; double b_n; double c_n; double D_plane;
    vector<double> crack_plane = {a_n, b_n, c_n, D_plane};

    //1
    subcomplex(unsigned int subcomplex_id_new) { // constructor 1, simple
        subcomplex_id = subcomplex_id_new;
    }
    //2
    subcomplex(unsigned int subcomplex_id_new, std::vector <unsigned int> new_sub_grains_sequence) { // constructor 2, based on a sub_grains_sequence
        subcomplex_id = subcomplex_id_new;
        Set_grains_sequence(new_sub_grains_sequence);
    }

    /// Grains
    void Set_grains_sequence(std::vector <unsigned int> new_sub_grains_sequence){
        sub_grains_sequence = new_sub_grains_sequence; }

    std::vector <unsigned int> Get_grains_sequence(unsigned int subcomplex_id){
        if(sub_grains_sequence.size() != 0)
            return sub_grains_sequence; }


    /// geometry
    void Set_sub_grain_coordinates(std::vector<tuple<double, double, double>> new_sub_grain_coordinates){
        sub_grain_coordinates = new_sub_grain_coordinates; }
    std::vector<tuple<double, double, double>> Get_sub_grain_coordinates(unsigned int subcomplex_id){
        return sub_grain_coordinates; }

    /* //2
    std::vector <unsigned int>  Get_grains_sequence(unsigned int subcomplex_id, std::vector <unsigned int> new_sub_grains_sequence){
        if(sub_grains_sequence.size() != 0) return sub_grains_sequence;
        else { Set_grains_sequence(new_sub_grains_sequence); return sub_grains_sequence; } }
        */

    /// Faces
    void Set_faces_sequence(std::vector <unsigned int> new_sub_faces_sequence){
        sub_faces_sequence = new_sub_faces_sequence; }
    //1
    std::vector <unsigned int>  Get_faces_sequence(unsigned int  subcomplex_id){
        if(sub_faces_sequence.size() != 0)
            return sub_faces_sequence; }

    void Set_common_faces_sequence(std::vector <unsigned int> new_common_faces_sequence){
        common_faces_sequence = new_common_faces_sequence; }
    std::vector <unsigned int> Get_common_faces_sequence(unsigned int subcomplex_id){
        return common_faces_sequence; }

/*    //2
    std::vector <unsigned int>  Get_faces_sequence(unsigned int  subcomplex_id, std::vector <unsigned int> new_sub_faces_sequence){
        if(sub_faces_sequence.size() != 0)
            return sub_faces_sequence;
        else {
            Set_grains_sequence(new_sub_faces_sequence);
            return sub_faces_sequence; }
    } */

    void Set_sfaces_sequence(std::vector <unsigned int> const &ssub_faces_sequence){
        s_sub_faces_sequence = ssub_faces_sequence;
    }
    std::vector <unsigned int> Get_sfaces_sequence(unsigned int  subcomplex_id){
        if(s_sub_faces_sequence.size() > 0) return s_sub_faces_sequence;
        else cout << "Caution! s_sub_faces_sequence = 0! Please add it to the corresponding subcomplex/crack"s << endl;
    }

    void Set_cfaces_sequence(std::vector <unsigned int> csub_faces_sequence){
        c_sub_faces_sequence = csub_faces_sequence;
    }
    std::vector <unsigned int> Get_cfaces_sequence(unsigned int  subcomplex_id){
        return c_sub_faces_sequence;
    }

    /// geometry
    void Set_common_faces_coordinates(std::vector<tuple<double, double, double>> new_common_faces_coordinates){
        common_faces_coordinates = new_common_faces_coordinates; }

    std::vector<tuple<double, double, double>> Get_common_faces_coordinates(unsigned int subcomplex_id){
        return common_faces_coordinates; }

    subcomplex Get_half_plane(subcomplex new_sub, double crack_length, std::vector<unsigned int> const &half_sub_sfaces_sequence){
//        vector<tuple<double, double, double>> grain_coordinates = grain_coordinates_vector;
//REPAIR        cout << " grain_coordinates_vector " << grain_coordinates.size() << endl;
        //cout << "Xmin " << get<0>(minmax_tuple.at(0)) << " Ymin " <<get<1>(minmax_tuple.at(0)) << " Zmin " << get<2>(minmax_tuple.at(0)) << endl;
        vector<tuple<double, double, double>> face_coordinates = face_coordinates_vector;
//REPAIR        cout << " face_coordinates_vector " << face_coordinates.size() << endl;

        subcomplex half_plane_cut(1);
        std::vector<unsigned int> half_sub_grains_sequence, half_common_faces_sequence;
//        std::vector<unsigned int> sub_grains_sequence = half_sub.Get_grains_sequence(0);
//TEMPORARILY        std::vector<unsigned int> sub_common_face_sequence = half_plane_cut.Get_common_faces_sequence(0);
//TEMPORARILY        std::vector<unsigned int> sub_sfaces_sequence = half_plane_cut.Get_sfaces_sequence(0);
//TEMPORARILY         std::vector<tuple<double, double, double>> sub_common_face_coordinates = half_plane_cut.Get_common_faces_coordinates(0);
/* ARCHIVE
        for (auto  itr = grain_coordinates.begin(); itr != grain_coordinates.end(); ++itr)
            if (std::find(sub_grains_sequence.begin(), sub_grains_sequence.end(), distance(grain_coordinates.begin(),itr)) != sub_grains_sequence.end()
                    && get<0>(*itr) < crack_length) half_sub_grains_sequence.push_back(distance(grain_coordinates.begin(),itr));
*/
//REPAIR cout << "half_sub.Get_grains_sequence(0) " << half_sub.Get_grains_sequence(0).size() << endl;
        for (auto half_grains : new_sub.Get_grains_sequence(0)) {
            if (get<1>(grain_coordinates_vector.at(half_grains)) < crack_length) {
//REPAIR                cout << " half_grains " << half_grains << " half_grain_coordinates " << get<0>(grain_coordinates_vector.at(half_grains)) << " crack_length " << crack_length << endl;
                half_sub_grains_sequence.push_back(half_grains); //// get<0> ->> only along X now!!!
            }
            }
        /* TEMPORARILY!!!
        for (auto  itr2 = face_coordinates.begin(); itr2 != face_coordinates.end(); ++itr2)
            if (std::find(sub_common_face_sequence.begin(), sub_common_face_sequence.end(), distance(face_coordinates.begin(),itr2)) != sub_common_face_sequence.end()
                    && get<0>(*itr2) < crack_length) half_common_faces_sequence.push_back(distance(face_coordinates.begin(),itr2));

        for (auto  itr3 = face_coordinates.begin(); itr3 != face_coordinates.end(); ++itr3)
            if (std::find(sub_sfaces_sequence.begin(), sub_sfaces_sequence.end(), distance(face_coordinates.begin(),itr3)) != sub_sfaces_sequence.end()
                && get<0>(*itr3) < crack_length) half_sub_sfaces_sequence.push_back(distance(face_coordinates.begin(),itr3));
*/
        half_plane_cut.Set_grains_sequence(half_sub_grains_sequence); // all grains before cut
 ///       half_plane_cut.Set_common_faces_sequence(half_common_faces_sequence); // common faces
        half_plane_cut.Set_sfaces_sequence(half_sub_sfaces_sequence); //special faces
        //half_plane_cut.Set_common_faces_coordinates(common_faces_coordinates);
        //half_plane_cut.Set_sub_grain_coordinates(subcomplex_grain_coordinates);
        //half_plane_cut.Set_cfaces_sequence(c_sub_faces_sequence); //cracked (induced) faces
        half_plane_cut.sub_length = crack_length;

        return half_plane_cut;
    }

/*    //int axis_id; // id of a considered axis: 1 for x, 2 for y and 3 for z
    switch (axis_id) {
        case 1: {
            double edge_coord = length_ratio*Lx_size; break; }
        case 2: {
            double edge_coord = length_ratio*Ly_size; break; }
        case 3: {
            double edge_coord = length_ratio*Lz_size; break; }
    }
*/

}; /// end of class SUBCOMPLEX


/// #6# The class of a MACROCRACK
class macrocrack {
    double total_fracture_energy = 0;
    subcomplex half_plane_subcomplex; // geometry part

private:
    double a_n;
    double b_n;
    double c_n;
    double D_plane;
    double crack_length;
    double real_crack_length;
    vector<double> crack_plane = {a_n, b_n, c_n, D_plane};

public:
    int crack_id = 0;
    double surface_energy = 0;
    double bridging_energy = 0;
    double multiple_cracking_energy = 0;
    double stress_concentrators_energy = 0;

    //3
    macrocrack(int crack_id_new, subcomplex &half_plane_sub) : half_plane_subcomplex(0) { // constructor 3, with a subcomplex
        crack_id = crack_id_new;
        crack_length = half_plane_sub.sub_length; // crack length from the corresponding subcomplex size
        half_plane_subcomplex = half_plane_sub; //set subcomplex
    }

    double Get_crack_length(int crack_id_new) {
        return crack_length;
    }

    void Set_real_crack_length(double sample_size) {
        real_crack_length = half_plane_subcomplex.sub_length * sample_size;
    }
    double Get_real_crack_length() {
        if (real_crack_length > 0.0) return real_crack_length;
        else cout << "Caution! real_crack_length = 0! Please use Set_real_crack_length(double sample_size) before!"s << endl;
    }

    void Set_crack_plane() {
        crack_plane = {half_plane_subcomplex.a_n, half_plane_subcomplex.b_n, half_plane_subcomplex.c_n, half_plane_subcomplex.D_plane};
    }
    vector<double> Get_crack_plane() {
        if (crack_plane.size() > 0.0) return crack_plane;
        else cout << "Caution! crack_plane.size() = 0! Please update first {half_plane_subcomplex.a_n, half_plane_subcomplex.b_n, half_plane_subcomplex.c_n, half_plane_subcomplex.D_plane} in the corresponding subcomplex!"s << endl;
    }

    void Set_multiple_cracking_energy(double total_energy) {
        multiple_cracking_energy = total_energy;
    }
    double Get_multiple_cracking_energy() {
        return multiple_cracking_energy;
    }

    std::vector <unsigned int> Get_faces_sequence(){
            return half_plane_subcomplex.Get_faces_sequence(crack_id); }

    std::vector <unsigned int> Get_sfaces_sequence(){
        return half_plane_subcomplex.Get_sfaces_sequence(0); }

    std::vector <tuple<double,double,double>> Get_common_faces_coordinates(unsigned int  crack_id){
        return half_plane_subcomplex.Get_common_faces_coordinates(crack_id); }
    /*
    double Get_surface_energy()
    {
        if (crack_length != 0.0)
            return ;
        else
            return 0.0;
    }
*/

}; /// end of class MACROCRACK

#endif
#endif

