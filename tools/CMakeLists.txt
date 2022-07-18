message(STATUS "Preparing lepto-utils..")

if (${HEPGEN_DISABLE_ROOT} MATCHES "YES")

	message("*** WANT DISABLED ROOT - No lumiCalc added ***")
	target_link_libraries(hepgen)

else()

	if (ROOT_FOUND)
		add_executable(makePlotsFromFile makePlotsFromFile.cc)
		add_executable(lumiCalc lumiCalc.cc)
		add_executable(beamFilePlotter beamfileplotter.cc)
		add_executable(leptoFileAnalyzer leptoFileAnalyzer.cpp)
	
		target_link_libraries(beamFilePlotter ${ROOT_LIBRARIES} hepgen)
		target_link_libraries(makePlotsFromFile ${ROOT_LIBRARIES})
		target_link_libraries(lumiCalc ${ROOT_LIBRARIES} hepgen)
		target_link_libraries(leptoFileAnalyzer slread hepgen ${ROOT_LIBRARIES})
	
		install (TARGETS beamFilePlotter DESTINATION bin)
		install (TARGETS lumiCalc DESTINATION bin)
		install (TARGETS leptoFileAnalyzer DESTINATION bin)
	endif()
endif()


if (${USE_EXPERIMENTAL} MATCHES TRUE)
	add_executable(xSecPlotter xSecPlotter.cc)
	target_link_libraries(xSecPlotter hepgen Moutarde VGG cln)
endif()

add_executable(beamFileSplitter beamFileSplitter.cpp)
target_link_libraries(beamFileSplitter hepgen)


add_library(slread lfread.cpp lfwrite.cpp)


add_executable(leptoFileSplitter splitter.cpp)
target_link_libraries(leptoFileSplitter slread)


install (TARGETS beamFileSplitter DESTINATION bin)
install (TARGETS slread DESTINATION lib)
install (TARGETS leptoFileSplitter DESTINATION bin)
