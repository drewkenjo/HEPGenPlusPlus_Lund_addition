message(STATUS "Preparing library libGKPi0.so..")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

FILE(GLOB gkpi0list *.hh)

add_library(GKPi0 SHARED libGKPi0.cpp gkSubProcessTable.cpp)
target_link_libraries(GKPi0 ${ROOT_LIBRARIES} MathMore ${GSL_LIBRARIES})

add_executable(testGKPi0 gkTest.cpp)
target_link_libraries(testGKPi0 GKPi0 ${ROOT_LIBRARIES})

add_executable(gkCrossSection gkCrossSection.cpp)
target_link_libraries(gkCrossSection GKPi0 ${ROOT_LIBRARIES})

add_executable(tableNator tableNator.cpp)

install (TARGETS GKPi0 DESTINATION lib)
install (TARGETS gkCrossSection DESTINATION bin)
install (FILES ${gkpi0list} DESTINATION include)
#install (FILES ${gkpi0listData} DESTINATION share/gkpi0data)
