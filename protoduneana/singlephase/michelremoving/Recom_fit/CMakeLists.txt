cet_make_exec( Recombination_fit
               SOURCE protoDUNE_Recom_fit.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               CaloTools
	       cetlib::cetlib 
	       cetlib_except::cetlib_except
	       fhiclcpp::fhiclcpp
               messagefacility::MF_MessageLogger
               )
