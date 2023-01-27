///================================ DCC Writer module ===================================================================================///
///=====================================================================================================================================///
/** The function in this library perform a structured output of the calculation results **/
///===================================================================================================================================///

/// Attached user defined C++ libraries:
///-------------------------------------
#include "Writer_Functions.h"
//#include "DCC_Writer/Writer_Functions.h"
///-------------------------------------

void DCC_Writer(std::vector<unsigned int> const  &sfaces_sequence, std::vector<unsigned int> const  &kface_sequence, macrocrack const &large_crack, vector <unsigned int> &State_sVector, double const &mu_f, double const &sigm_f, vector<vector<int>> const &RW_series_vector,  const string P_type) {

    int output_counter = 0; // special counter for output numeration
    vector<agglomeration> agglomerations_vector; // vector containing aglomarations as its objects - only for L type simulations (!)


    char* stype = const_cast<char*>(P_type.c_str()); // Processing type variable
    char P_simulation_type = *stype;

    /// Output special and ordinary face sequences to the output directory specified in config.txt
    DCC_sequences_Writer(sfaces_sequence, kface_sequence, output_counter);
    if (P_simulation_type == 'L') {
        Strips_distribution_Writer(output_counter); // strips distribution
        RWseries_Writer(RW_series_vector, output_counter); // vector of all strips

        /// agglomerations  output
        Agglomerations_Writer(output_counter, RW_series_vector, mu_f, sigm_f);

    } // end of if(P_simulation_type == 'L')
     return ;
}