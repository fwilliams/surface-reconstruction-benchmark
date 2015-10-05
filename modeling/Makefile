TEEM_FLAG=0

CXX=g++

OPENMESH_INC = ../OpenMesh
ANN_INC = ../ann/include

OPENMESH_LIB = ../OpenMesh/lib
ANN_LIB = ../ann/lib

CXXFLAGS = -O3 -I$(OPENMESH_INC) -I$(ANN_INC) -g
LDFLAGS = -O3 -lblas -llapack -L$(OPENMESH_LIB) -lOpenMeshCore -L$(ANN_LIB) -lANN

CXXFLAGS += -DHAS_TEEM=0

ifeq ($(TERM_PROGRAM), Apple_Terminal)
	CXXFLAGS += -DHAS_MALLOC=0
else
	CXXFLAGS += -DHAS_MALLOC=1
endif

LIB_DIR=lib
BIN_DIR=../bin

LINALG_OBJS = $(LIB_DIR)/DenseMatrix.o $(LIB_DIR)/Vector.o $(LIB_DIR)/DenseEigenSystem.o $(LIB_DIR)/lm_method.o \
				  $(LIB_DIR)/RotationMatrix.o $(LIB_DIR)/Rigid.o
IMPLICIT_OBJS = $(LIB_DIR)/ImplicitFunction.o $(LIB_DIR)/mpu_soup.o $(LIB_DIR)/tri_integrator.o

LIB_OBJS = $(LINALG_OBJS) $(IMPLICIT_OBJS)
MESH_IMP_OBJ = $(LIB_DIR)/mesh_to_implicit.o
ISOSURFACE_OBJ = $(LIB_DIR)/isosurface.o

default : all
all : lib mesh_imp isosurface

lib : $(LIB_OBJS)
	@ echo "Creating static library for modeling ..."
	@ ar cr $(LIB_DIR)/libmodeling.a $(LIB_OBJS)
		@ ranlib $(LIB_DIR)/libmodeling.a

mesh_imp : $(LIB_OBJS) $(MESH_IMP_OBJ)
	@ echo "creating mesh to implicit..."
	@ $(CXX) -o $(BIN_DIR)/mesh_to_implicit $(LIB_OBJS) $(MESH_IMP_OBJ) $(LDFLAGS)

isosurface : $(LIB_OBJS) $(ISOSURFACE_OBJ)
	@ echo "creating isosurfacer..."
	@ $(CXX) -o $(BIN_DIR)/isosurface $(LIB_OBJS) $(ISOSURFACE_OBJ) $(LDFLAGS) 

$(LIB_DIR)/%.o : %.cpp
	@ echo "Compiling $< ..."
	@ $(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	rm -f $(LIB_DIR)/*.o $(LIB_DIR)/libmodeling.a $(BIN_DIR)/mesh_to_implicit $(BIN_DIR)/isosurface
