#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=MinGW_2014-Windows
CND_DLIB_EXT=dll
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/E_LJ.o \
	${OBJECTDIR}/E_REBO_C.o \
	${OBJECTDIR}/E_REBO_CH.o \
	${OBJECTDIR}/E_Torsion.o \
	${OBJECTDIR}/OpenCL_tools.o \
	${OBJECTDIR}/OpenREBO.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/nlist.o \
	${OBJECTDIR}/splines.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-ffast-math -DNDEBUG -march=native -pg
CXXFLAGS=-ffast-math -DNDEBUG -march=native -pg

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/openrebo.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/openrebo.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/openrebo ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/E_LJ.o: E_LJ.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/E_LJ.o E_LJ.cpp

${OBJECTDIR}/E_REBO_C.o: E_REBO_C.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/E_REBO_C.o E_REBO_C.cpp

${OBJECTDIR}/E_REBO_CH.o: E_REBO_CH.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/E_REBO_CH.o E_REBO_CH.cpp

${OBJECTDIR}/E_Torsion.o: E_Torsion.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/E_Torsion.o E_Torsion.cpp

${OBJECTDIR}/OpenCL_tools.o: OpenCL_tools.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OpenCL_tools.o OpenCL_tools.cpp

${OBJECTDIR}/OpenREBO.o: OpenREBO.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OpenREBO.o OpenREBO.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/nlist.o: nlist.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/nlist.o nlist.cpp

${OBJECTDIR}/splines.o: splines.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/splines.o splines.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/openrebo.exe

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
