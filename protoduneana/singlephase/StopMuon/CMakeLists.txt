art_make( 
  LIBRARY_NAME ProtoDUNEStopMuonLib 
  LIB_LIBRARIES 
  ProtoDUNEDataUtils
  ProtoDUNEUtilities
  dune_DuneObj
  lardataobj_RecoBase
  lardataalg_detectorinfo

	  larcore_Geometry_Geometry_service
	  lardata_Utilities
	  larevt_Filters
	  larreco_RecoAlg

  MODULE_LIBRARIES
  larcorealg_Geometry
	lardataobj_RawData
	lardataobj_RecoBase
  lardataobj_AnalysisBase #EXTRA?
	lardata_Utilities
	larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
  larsim_MCCheater_BackTrackerService_service #EXTRA
  larsim_MCCheater_ParticleInventoryService_service #EXTRA
  lardata_RecoObjects #EXTRA
  larreco_RecoAlg #EXTRA
  larevt_Filters #EXTRA
  larpandora_LArPandoraInterface #EXTRA

  ProtoDUNEStopMuonLib 
  ProtoDUNEDataUtils
  ProtoDUNEUtilities
  ProtoDUNEStopMuonLib 

  lardataalg_detectorinfo
   
	larcore_Geometry_Geometry_service

  lardata_ArtDataHelper
  lardataalg_detectorinfo
	SignalShapingServiceDUNE10kt_service
	SignalShapingServiceDUNE35t_service
  nusimdata_SimulationBase
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
	
	${FHICLCPP}
	cetlib cetlib_except
  ${ROOT_GEOM}
  ${ROOT_XMLIO}
	${ROOT_BASIC_LIB_LIST}

  )

install_headers()
install_fhicl()
install_source()
#install_fw( LIST toyWaveform.txt)
