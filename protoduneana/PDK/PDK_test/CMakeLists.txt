art_make(BASENAME_ONLY
	 	
	 MODULE_LIBRARIES
	 #protoduneana/Utilities
         ProtoDUNEUtilities
         lardataalg_DetectorInfo
	 lardataobj_RawData
	 lardataobj_RecoBase
	 larcorealg_Geometry
	 larcore_Geometry_Geometry_service
	 lardata_Utilities
	 larsim_MCCheater_BackTrackerService_service
	 larsim_MCCheater_ParticleInventoryService_service
	 nusimdata_SimulationBase
         lardata_ArtDataHelper
         lardataobj_Simulation
         larreco_Calorimetry
         larevt_Filters
         larreco_RecoAlg
         larpandora_LArPandoraInterface
       	 ${LARRECO_LIB}
	 ${LARDATA_LIB}
	 ${ART_FRAMEWORK_CORE}
	 ${ART_FRAMEWORK_PRINCIPAL}
	 ${ART_FRAMEWORK_SERVICES_REGISTRY}
	 ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
	 ${ART_ROOT_IO_TFILESERVICE_SERVICE}
	 art_Persistency_Common canvas
	 art_Persistency_Provenance canvas
	 art_Utilities canvas
	 ${MF_MESSAGELOGGER}
	 
	 cetlib cetlib_except 
	 LIB_LIBRARIES
         ${ROOT_BASIC_LIB_LIST}
	 ${ROOT_GEOM}
	 ${ROOT_XMLIO}
         ${ROOT_GDML}
)
art_make_exec(NAME Template_PIDloglike_range
SOURCE Template_PIDloglike_range.cpp
LIBRARIES
protoduneana_PDK
${ART_FRAMEWORK_CORE}
${ART_FRAMEWORK_PRINCIPAL}
${ART_FRAMEWORK_SERVICES_REGISTRY}
${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
${ART_ROOT_IO_TFILESERVICE_SERVICE}
${ROOT_BASIC_LIB_LIST}
${ROOT_GEOM}
${ROOT_XMLIO}
${ROOT_GDML} 
)


#add_subdirectory(include)
#add_subdirectory(BDT)

install_headers()
install_fhicl()
install_scripts()
