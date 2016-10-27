# ----------------------
# Configuration options
# ----------------------
INCLUDE(CMakeDependentOption)

OPTION(ENABLE_HOMEWORK "Enable Chrono ME751 Homework Problems" ON)

IF(NOT ENABLE_HOMEWORK)
    RETURN()
ENDIF()

#--------------------------------------------------------------
# List here the names of all tests

MESSAGE(STATUS "Adding Homework Problems...")

SET(HW_PROGRAMS
    ME751_HW6_P4
	ME751_HW7_P2
)

#--------------------------------------------------------------

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  SET(WORK_DIR ${PROJECT_BINARY_DIR}/bin/$<CONFIGURATION>)
ELSE()
  SET(WORK_DIR ${PROJECT_BINARY_DIR}/bin)
ENDIF()

#--------------------------------------------------------------
# Add executables

FOREACH(PROGRAM ${HW_PROGRAMS})
  MESSAGE(STATUS "... ${PROGRAM}")
  
  ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
  SOURCE_GROUP(""  FILES  "${PROGRAM}.cpp")

  SET_TARGET_PROPERTIES(${PROGRAM}  PROPERTIES
    FOLDER tests
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    LINK_FLAGS "${COMPILE_FLAGS} ${CHRONO_LINKER_FLAGS}"
    )

  TARGET_LINK_LIBRARIES(${PROGRAM} ${CHRONO_LIBRARIES})

  # Note: this is not intended to work on Windows!
  ADD_TEST(NAME ${PROGRAM}
           WORKING_DIRECTORY ${WORK_DIR}
           COMMAND ${WORK_DIR}/${PROGRAM}
           )

ENDFOREACH()

