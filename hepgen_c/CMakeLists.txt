message(STATUS "Preparing HEPGen++ binary")

set(linklist hepgen)
set(CURSES_NEED_NCURSES TRUE)

if (ROOT_FOUND)
  set(linklist ${linklist} ${ROOT_LIBRARIES})
  set(USE_ROOT 1)
else()
  message(STATUS "NO ROOT FOUND! Building ultra-basic version of libhepgen.so")
endif()

if (CURSES_FOUND)
  add_executable(HEPGenCurses++ main_curses.cpp)
  target_link_libraries(HEPGenCurses++ ${linklist} ${CURSES_LIBRARIES})
  install (TARGETS HEPGenCurses++ RUNTIME DESTINATION bin)
endif()

link_directories(${CMAKE_BINARY_DIR})

add_executable(HEPGen++ main.cpp)
target_link_libraries(HEPGen++ ${linklist})

install (TARGETS HEPGen++ RUNTIME DESTINATION bin)
