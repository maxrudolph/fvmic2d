#Makefile for FEM code
MARKERCODE_HG_VERSION = `hg tip | grep changeset | cut -d : -f 2,3 | tr -d \ `
#PETSC_DIR	= /usr/local/petsc-2.3.3-p13
#PETSC_ARCH	= linux-gnu-c-debug
#CC = gcc-4.2
PETSC_INCLUDE	= #${PETSC_DIR}/include
FEM_INCLUDE =   #../include ../src ../src/Elements
#VPATH	=	../src ../include ../src/Elements ../tests ../src/Materials ../src/Thermal ../src/CoupledThermalDamage ../src/Damage/
#LOCDIR		= /home/max/projects/fem 
#LOCDIR		= /Users/max/projects/fem
#LOCDIR		= /Users/max/projects/markercode/vep/ridges/parallel/
LOCDIR          = /data/64/max/projects/markercode/vep/ridges/parallel/
#CFLAGS 	= -Wall -O3  -I${PETSC_DIR}/src/dm/mesh/sieve ${BOOST_INCLUDE} ${TRIANGLE_INCLUDE} ${TETGEN_INCLUDE} ${PETSC_INDLUDE} ${FEM_INCLUDE}
CFLAGS		= -Wall -O0 -g -fsanitize=undefined -fno-omit-frame-pointer
#CFLAGS 	= -Wall -O3  -I${PETSC_DIR}/src/dm/mesh/sieve ${BOOST_INCLUDE} ${TRIANGLE_INCLUDE} ${TETGEN_INCLUDE} ${PETSC_INDLUDE} ${FEM_INCLUDE}
#CFLAGS 	= -Wall -O0  -I${PETSC_DIR}/src/dm/mesh/sieve ${BOOST_INCLUDE} ${TRIANGLE_INCLUDE} ${TETGEN_INCLUDE} ${PETSC_INDLUDE} ${FEM_INCLUDE}
FFLAGS	=
#LIBDIR  = lib
CPPFLAGS 	=
#
NP		= 1
MANSEC		= KSP

#include ${PETSC_DIR}/bmake/common/base
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

all: subduction

version.h : 
	echo "#define MARKERCODE_HG_VERSION \"$(MARKERCODE_HG_VERSION)\"" > version.h

fdcode.h : version.h

io.o : version.h

clean ::
	rm -f *.o version.h subduction



#Fem : main.o linearHexahedralElement.o matrixops.o gaussRules.o input.o mesh.o
#Fem : main.o taylorHoodElement.o matrixops.o gaussRules.o input.o mesh.o elementMassArray.o nonLinearElastic.o
#fddriver : main.o kspLinearSolve.o
#	-${CLINKER} $(CFLAGS) -o fddriver $^ ${PETSC_LIB}
#	${RM} main.o

#vep-ridge : main_vep_ridges.o kspLinearSolve.o materialProperties.o markers.o nodalFields.o formThermalSystem.o reggrid.o nodalStressStrain.o vep_system.o updateMarkerStrainPressure.o updateMarkerStress.o subgridStressChanges.o io.o updateDamageViscosity.o options.o gridMarkers.o
#	-${CLINKER} $(CFLAGS) -o vep-ridge $^ ${PETSC_LIB}
#	${RM} main.o

main_subduction.o : main_subduction.c version.h
	$(CC) $(CFLAGS) -c $^ $(PETSC_CC_INCLUDES)


subduction:  post.o main_subduction.o nodalFields.o gridGenerator.o gridSpacing.o options.o  markers.o randomProperties.o io.o kspLinearSolve.o vep_system.o retrieveSolutions.o limitTimestep.o nodalStressStrain.o updateMarkerStrainPressure.o subgridStressChanges.o updateMarkerStress.o formThermalSystem.o subgridTemperatureChanges.o updateBCs.o adiabaticHeating.o enforceThermalBC.o restart.o viscosity.o benchmarkInitialConditions.o residual.o profile.o markerProjection.o initialPressureGuess.o pressureNullSpace.o
	-${CLINKER} $(CFLAGS) -o subduction $^ ${PETSC_LIB}
	${RM} main.o


#To compile with texture enabled, run:
#make clean; make CFLAGS="-DTEXTURE" testMM
testMM : testMicroMacro.o microMacroIceModel.o tensorOps.o
	-${CLINKER} $(CFLAGS) -o testMM $^ ${PETSC_LIB}
	${RM} main.o

initializeFabric : initializeFabricFromMarkers.o
	-${CLINKER} $(CFLAGS) -o initializeFabric $^ ${PETSC_LIB}
	${RM} main.o

consolidateTexture : consolidateTexture.o
	-${CLINKER} $(CFLAGS) -o consolidateTexture $^ ${PETSC_LIB}
	${RM} main.o

testls : main_testlogsummary.o
	-${CLINKER} $(CFLAGS) -o testls $^ ${PETSC_LIB}
	${RM} main.o



comparemmgk : comparemmgk.o microMacroIceModel.o tensorOps.o viscosity.o
	-${CLINKER} $(CFLAGS) -o testMM $^ ${PETSC_LIB}
	${RM} main.o

gridMarkers : gridMarkersPost.o
	-${CLINKER} $(CFLAGS) -o gridMarkers $^ ${PETSC_LIB}
	${RM} main.o

textureIndex : texture_index.o
	-${CLINKER} $(CFLAGS) -o textureIndex $^ ${PETSC_LIB}
	${RM} main.o



testviscosity : test_viscosity.o viscosity.o
	-${CLINKER} $(CFLAGS) -o testviscosity $^ ${PETSC_LIB}
	${RM} main.o

testPhaseBA : testPhaseBA.o phaseBA.o 
	-${CLINKER} $(CFLAGS) -o testPhaseBA $^ ${PETSC_LIB}
	${RM} main.o

#testgps : testGps.o
#	-${CLINKER} $(CFLAGS) -o testgps $^ ${PETSC_LIB}
#	${RM} main.o
