# art_make assumes all non-plugin source files go into a library.
#
# Dictionary is also made since classes_def.xml exists.
art_make(
  LIB_LIBRARIES
  bprinter
  ${ARTDAQ_DAQDATA} 
  ${ROOT_HIST}
	${ROOT_TREE}
	   ${ROOT_RINT}
   ${ROOT_THREAD}
   ${ROOT_CORE}
   ${ROOT_MATHCORE}
   ${ROOT_MATRIX}
   ${DS50DAQ_ONLINE_ARTSERVICES}
	NO_INSTALL
  )



install_headers(SUBDIRS detail)
install_source()
