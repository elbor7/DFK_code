///================================ A part of the DCC Multiphysics module =============================================================///
///=================================================================================================================================///
/** The library contains functions saetting various physical quantities for k-Cells of PCC               *
 *                                                                                              **/
///================================================================================================================================///
using namespace std;

/// Just three overloaded templates for "sign()" function [taken from the internet]
template <typename TP> inline constexpr
int sign(TP x, std::false_type is_signed) {
    return TP(0) < x;
}

template <typename TP> inline constexpr
int sign(TP x, std::true_type is_signed) {
    return (TP(0) < x) - (x < TP(0));
}

template <typename TP> inline constexpr
int sign(TP x) {
    return sign(x, std::is_signed<TP>());
}
//------------------------------------------------------------------------///
void crack_modes_stress_field(std::vector <double> &face_energies, int crack_mode, double new_crack_length, double external_vonMizes_stress, double Puasson_coeff = 0.3, double Young_modulus = 215.0*pow(10.0,9)) {
    double nu = Puasson_coeff;
    double GB_width = 3.0*pow(10.0,-9);
    double Len = new_crack_length;
    double Sigm = external_vonMizes_stress;
    double Sxx = 0.0, Syy = 0.0, Szz = 0.0, Sxy = 0.0;

    for (unsigned int itr = 0; itr < face_coordinates_vector.size()-1; ++itr) { // loop over all GBs and their barycentic coordinates
        double z = get<0>(face_coordinates_vector.at(itr))*sample_size_vector.at(0);
        double x = get<1>(face_coordinates_vector.at(itr))*sample_size_vector.at(1);
        double y = get<2>(face_coordinates_vector.at(itr))*sample_size_vector.at(2);

        switch (crack_mode) {
            case 1: { // crack mode I
                /// The first crack mode (Syy only!) ::
                double r2p = pow(x, 2.0) + pow(y, 2.0);
                double r2m = pow(x, 2.0) - pow(y, 2.0);
                double x2p = pow(x, 2.0) + pow(Len, 2.0);
                double x2m = pow(x, 2.0) - pow(Len, 2.0);
                double y2p = pow(x, 2.0) + pow(Len, 2.0);
                double y2m = pow(x, 2.0) - pow(Len, 2.0);
                double Pp = std::sqrt(pow(r2m, 2.0) + 4.0*pow(x, 2.0)*pow((y - Len),2));
                double Qp = std::sqrt(Pp + r2m);
                double Qm = std::sqrt(Pp - r2m);
                Sxx = Sigm * (Pp*r2p*(abs(y) * Qm - abs(x) * Qp) - abs(x) * Qp*( pow(y, 4.0) - 2.0*x2p*pow(y, 2.0) - 3.0*pow(x2m, 2.0) ) + abs(y) * Qm*( pow(y, 4.0) + 6.0*x2p*pow(y, 2.0) + 5.0*pow(x2m, 2.0) ) ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
                Syy = Sigm * (Pp*r2p*(abs(x) * Qp - abs(y) * Qm) + abs(x) * Qp*( 5.0*pow(y, 4.0) + 6.0*x2p*pow(y, 2.0) + pow(x2m, 2.0) ) + abs(y) * Qm*( 3.0*pow(y, 4.0) + 2.0*x2p*pow(y, 2.0) - pow(x2m, 2.0) ) ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
                Sxy = Sigm * sign(x*y) * (Pp*r2p*(abs(x) * Qm + abs(y) * Qp) - abs(x) * Qm*( 3.0*pow(y, 4.0) + 2.0*x2p*pow(y, 2.0) - pow(x2m, 2.0) ) + abs(y) * Qp*( pow(y, 4.0) - 2.0*x2p*pow(y, 2.0) - 3.0*pow(x2m, 2.0) ) ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
/*
                double Pp = std::sqrt(pow(x, 2.0) - pow(y, 2.0) - pow(pow(Len, 2) / 4.0, 2.0) + 4.0 * pow(x, 2.0) * pow(y, 2.0));
                double Qp = std::sqrt(Pp + (pow(x, 2.0) - pow(y, 2.0) - pow(Len, 2.0) / 4.0));
                double Qm = std::sqrt(Pp - (pow(x, 2.0) - pow(y, 2.0) - pow(Len, 2.0) / 4.0));

                Sxx = Sigm * (
                        Pp * (pow(x, 2.0) + pow(y, 2.0) - pow(Len, 2.0) / 4.0) * (abs(y) * Qm - abs(x) * Qp) -
                        abs(x) * Qp * (pow(y, 4.0) - 2.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) -
                                       3.0 * pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0)) +
                        abs(y) * Qm * (pow(y, 4.0) + 6.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) +
                                       5.0 * pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0))
                ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
                Syy = Sigm * (
                        Pp * (pow(x, 2.0) + pow(y, 2.0) - pow(Len, 2.0) / 4.0) * (abs(x) * Qp - abs(y) * Qm) +
                        abs(x) * Qp * (5.0 * pow(y, 4.0) + 6.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) +
                                       pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0)) +
                        abs(y) * Qm * (3.0 * pow(y, 4.0) + 2.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) -
                                       pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0))
                ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
                //+ Sigm;
                Szz = nu * (Sxx + Syy);
                Sxy = Sigm * sign(x * y) *
                             (Pp * (pow(x, 2.0) + pow(y, 2.0) - pow(Len, 2.0) / 4.0) * (abs(x) * Qm + abs(y) * Qp) -
                              abs(x) * Qm *
                              (3.0 * pow(y, 4.0) + 2.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) -
                               pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0)) +
                              abs(y) * Qp * (pow(y, 4.0) - 2.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) -
                                             3.0 * pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0))
                             ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
                cout << " Young_modulus: " << pow(x, 2.0) - pow(y, 2.0) - pow(pow(Len, 2) / 4.0, 2.0) + 4.0 * pow(x, 2.0) * pow(y, 2.0) << " face_energies.at(itr): " << face_energies.at(itr) << endl;
*/
            }
        case 2: { // crack mode II
    /// The second and the third (Sxy only) crack modes ::
       // S2xx := 0.5*Sigm* Sign[x*y]*(Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[x]*Qm + Abs[y]*Qp) + Abs[x]*Qm*(1.0*y^4 + 6.0*(x^2 + Len^2/4.0)*y^2 + 5.0*(x^2 - Len^2/4.0)^2) - Abs[y]*Qp*(3.0*y^4 + 10.0*(x^2 + Len^2/4.0)*y^2 + 7.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
        //S2yy := -0.5*Sigm* Sign[x*y]*(Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[x]*Qm + Abs[y]*Qp) - Abs[x]*Qm*(3.0*y^4 + 2.0*(x^2 + Len^2/4.0)*y^2 - 1.0*(x^2 - Len^2/4.0)^2) + Abs[y]*Qp*(1.0*y^4 - 2.0*(x^2 + Len^2/4.0)*y^2 - 3.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
        //S2xy := -0.5* Sigm* (Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[y]*Qm - Abs[x]*Qp) - Abs[x]*Qp*(1.0*y^4 - 2.0*(x^2 + Len^2/4.0)*y^2 - 3.0*(x^2 - Len^2/4.0)^2) + Abs[y]*Qm*(1.0*y^4 + 6.0*(x^2 + Len^2/4.0)*y^2 + 5.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
         }
         case 3: { // crack mode III
    /// The third crack mode ::
        //dl = 0.0
        //xp := x + dl
        //xm = x - dl
        //Pp3 := ((xp^2 - y^2 + y*Len)^2 + 4.0*xp^2*(y - Len/2.0)^2)^0.5
        //Qm3 := (Pp3^2 - (xp^2 - y^2 + y*Len))^0.5
        //Qp3 := (Pp3^2 + (xp^2 - y^2 + y*Len))^0.5
        //S3xy := 0.5* Sigm* (Pp3^2*(xp^2 + y^2 - y*Len)*(Qp3*(xp - y + Len/2.0) + Qm3*(xp + y - Len/2.0)) + Qm3*(xp + y - Len/2.0)*(1.0*xp^4 - 2.0*xm^2*(y^2 - y*Len + Len^2/2.0) - 3.0*y^2*(y - Len)^2) + Qp3*(y - Len/2.0)*(3.0*xp^4 + 2.0*xm^2*(y^2 - y*Len + Len^2/2.0) - y^2*(y - Len)^2) + Qp3*xp*(xp^4 + 6.0*xm^2*(y^2 - y*Len + Len^2/2.0) + 5.0*y^2*(y - Len)^2))/(2.0*(2.0^(0.5))*Pp3^6)
                }
        } // end of switch

        double SigmMAX = (1.0 + nu) * (Sxx + Syy + Szz) / 3.0;

        /// von Mizes stress
        double res_vonMizes_stress = 0.0;
        if (x <= Len)  res_vonMizes_stress = SigmMAX; /// assumption for the maximal stress value (!)
            else res_vonMizes_stress = std::sqrt(0.5 * (pow((Sxx - Syy),2.0) + pow((Sxx - Szz),2.0) + pow((Syy - Szz),2.0) + 6 * pow(Sxy,2.0)));

        /// von Mizes energies
        /// (1) sample_size_vector.at(1) now!
        if (res_vonMizes_stress >  0.0) face_energies.at(itr) = pow(res_vonMizes_stress,2)*face_areas_vector.at(itr)*pow(sample_size_vector.at(1),2.0)*GB_width/Young_modulus;
//REPAIR        cout << " itr: " << itr << " face_energies.at(itr): " << face_energies.at(itr) << endl;

    } // end of for( itr )

} /// end of crack_modes_stress_field() function



/// HEAP
/*
std::vector<double> Energy_GBs_state() {
    std::vector<double> total_GBenergies_vector; // main function output
    std::vector<grain_boundary> all_GBs; // vector of all GBs as structures

    for (unsigned int m = 0; m < CellNumbs.at(2); m++)
        all_GBs.push_back(grain_boundary(m));

    face_areas_vector

    GB_SE_vector.at(m) = 2.0*_surface_energy*face_areas_vector.at(m);
    GB_EEE_vector.at(m) =
    //von Mizes stress
    double s_mis = sqrt(0.5*((s11 - s22)^2 + (s22 - s33)^2 + (s33 - s11)^2 + 6.0*(s12^2 + s23^2 + s31^2)));
    double ext_elastic_energy = pow(s_mis,2)*face_areas_vector.at(m)*hm/Young_modulus;

    GB_CIE_vector =
    GB_BLE_vector =
    GB_CLE_vector =

    for(unsigned int m = 0; m < CellNumbs.at(2); m++) {
        all_GBs.at(m).Set_surface_energy(GB_SE_vector);
        all_GBs.at(m).Set_external_elastic_energy(GB_EEE_vector);
        all_GBs.at(m).Set_crack_interaction_energy(GB_CIE_vector);
        all_GBs.at(m).Set_Bl_energy(GB_BLE_vector);
        all_GBs.at(m).Set_Cl_energy(GB_CLE_vector);
    } // end of for (unsigned int m = 0; m < CellNumbs.at(2); m++)

    for(unsigned int m = 0; m < CellNumbs.at(2); m++)
        total_GBenergies_vector.push_back(all_GBs.at(m).Get_surface_energy(m) + all_GBs.at(m).Get_external_elastic_energy(m) +
                                                  all_GBs.at(m).Get_crack_interaction_energy(m) + all_GBs.at(m).Get_Bl_energy(m) + all_GBs.at(m).Get_Bl_energy(m));
        return total_GBenergies_vector;
} // end of vector<double> Energy_GBs_state()

/*
void data_multiphysics_reader() {
    std::string line;
    ifstream in_data(config);

    if (in_data) { //If the file was successfully open, then
        while(getline(in_data, line, '\n'))
            for (auto it: line) {

                if (it == '?') {
            stringstream line1_stream(line);
            line1_stream >> p_max;
    } else cout << "SIMULATION_MODE() error: The file " << config << " cannot be read" << endl; // If something goes wrong


        // Number of types
        for (auto it: line) {
            if (it == '&') Processing_type = line.at(0); // simulation type
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


        }

        */