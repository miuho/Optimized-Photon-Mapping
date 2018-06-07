add_executable(p4 main.cpp raytracer.cpp)
target_link_libraries(p4 application math scene tinyxml ${SDL_LIBRARY}
                      ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES}
                      ${GLEW_LIBRARIES})

if(APPLE)
    target_link_libraries(p4 SDLmain)
endif()

install(TARGETS p4 DESTINATION ${PROJECT_SOURCE_DIR}/..)
