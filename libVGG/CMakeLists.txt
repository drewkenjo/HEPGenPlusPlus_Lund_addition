message(STATUS "Preparing library libVGG.so..")

FILE(GLOB vgglist *.hh)
FILE(GLOB vgglistData data/*)

add_library(VGG SHARED GPDQ.cc mrst99.cc ../libhepgen/myTHEO.cc)
target_link_libraries(VGG ${ROOT_LIBRARIES})

install (TARGETS VGG DESTINATION lib)
install (FILES ${vgglist} DESTINATION include)
install (FILES ${vgglistData} DESTINATION share/vggdata)
