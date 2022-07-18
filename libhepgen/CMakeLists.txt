message(STATUS "Preparing library libhepgen.so..")

####### make the headerlist for installation
FILE(GLOB headerlist *.h)
FILE(GLOB headerlist_tools tools/*.h)

#### add the mytheo hh stuff
set(headerlist
        ${headerlist}
        myTHEO.hh
        )

#### just the normal libs
set(liblist
    hvector.cpp
    hlorentzvector.cpp
    hparammanager.cpp
    hcardparser.cpp
    hhelper.cpp
    hpionicdata.cpp
    hbooker.cpp
    hbookbackendASCII.cpp
    hphysicsgen.cpp
    hgenmanager.cpp
    hevent.cpp
    hparticle.cpp
    hdvcsgen.cpp
    houtput.cpp
    houtputbackend.cpp
    houtputbackendLEPTO.cpp
    hrhogen.cpp
    hpigen.cpp
    houtputbackendLEPTOGF.cpp
    hbeamfile.cpp
    hphigen.cpp
    hrhoplusgen.cpp
    homegagenpigamma.cpp
    hjpsigen.cpp
    hrotmat.cpp
    homegagen.cpp
    hpamgen.cpp
    reweightKine.cc
    myTHEO.cc
    )


  set(liblist ${liblist}
    hbookbackendROOT.cpp
    houtputbackendROOT.cpp
    hMosseGen.cpp)
  set(linklist ${linklist} ${ROOT_LIBRARIES})
  set(USE_ROOT 1)

#### BUILD THE LIBRARY
add_library(hepgen SHARED ${liblist})

target_link_libraries(hepgen ${Geant4_LIBRARIES})

#### LINK TO ROOT
target_link_libraries(hepgen VGG)

#### INSTALLATION #####
install (FILES ${headerlist} DESTINATION include)
install (FILES ${headerlist_tools} DESTINATION include)
install (DIRECTORY CLHEP_EMBEDDED DESTINATION include)
install (TARGETS hepgen LIBRARY DESTINATION lib)
