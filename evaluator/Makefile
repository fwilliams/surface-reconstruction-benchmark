OPENMESH_INC = ../OpenMesh
OPENMESH_LIB = ../OpenMesh/lib
ANN_INC = ../ann/include

MODELING_LIB = ../modeling/lib
SAMPLER_LIB = ../sampler/lib

PBRT_INCS = -I../pbrt/src -I../pbrt/src/core -I../pbrt/src/accelerators -I../pbrt/src/materials -I../pbrt/src/shapes
PBRT_LIB = ../pbrt/src/objs
ANN_LIB = ../ann/lib

DEFS = -DNDEBUG -DPBRT_PROBES_NONE

CXX=g++
CXXFLAGS = -O3 -I$(OPENMESH_INC) $(PBRT_INCS) -I$(ANN_INC)
LDFLAGS = -O3 -L$(OPENMESH_LIB) -L$(MODELING_LIB) -L$(SAMPLER_LIB) -L$(PBRT_LIB) -L$(ANN_LIB)
CXX_LIBS = -lpng14 -lblas -llapack -lmodeling -lsampler -lOpenMeshCore -lpbrt -ldl -lANN

ifeq ($(TERM_PROGRAM), Apple_Terminal)
	CXXFLAGS += -DHAS_MALLOC=0
	CXX_LIBS += -lz
else
	CXXFLAGS += -DHAS_MALLOC=1
	CXX_LIBS += -lpthread
endif

LIB_DIR=lib
BIN_DIR=../bin

TEEM_FLAG=0
CXXFLAGS += -DHAS_TEEM=0 $(DEFS)

VPATH=.:./plot

EVAL_OBJS = $(LIB_DIR)/KdOpenMesh.o $(LIB_DIR)/kd_tree.o $(LIB_DIR)/mesh_sampler.o $(LIB_DIR)/ShortestDistanceMap.o \
				$(LIB_DIR)/recon_evaluation.o $(LIB_DIR)/shape_distribution.o $(LIB_DIR)/GlobalStats.o \
				$(LIB_DIR)/implicit_sampler.o $(LIB_DIR)/UniformGrid.o
LIB_OBJS = $(EVAL_OBJS)

RUNNER_OBJS = $(LIB_OBJS) $(LIB_DIR)/run_evaluation.o
SAMPLEMESH_OBJS = $(LIB_OBJS) $(LIB_DIR)/sample_mesh.o
FIXER_OBJS = $(LIB_OBJS) $(LIB_DIR)/mesh_fixer.o
SAMPLER_OBJS = $(LIB_OBJS) $(LIB_DIR)/implicit_uniform.o

SINGLE_OBJS = $(LIB_OBJS) $(LIB_DIR)/single_distribution.o
AGGREGATE_OBJS = $(LIB_OBJS) $(LIB_DIR)/aggregate_distribution.o

default : all
all : lib run sample_mesh implicit_uniform single aggregate

run : $(RUNNER_OBJS)
	@ echo "creating evaluation runner..."
	@ $(CXX) -o $(BIN_DIR)/run_evaluation $(RUNNER_OBJS) $(LDFLAGS) $(CXX_LIBS)

sample_mesh : $(SAMPLEMESH_OBJS)
	@ echo "creating mesh sampler..."
	@ $(CXX) -o $(BIN_DIR)/sample_mesh $(SAMPLEMESH_OBJS) $(LDFLAGS) $(CXX_LIBS)

fixer : $(FIXER_OBJS)
	@ echo "creating mesh_fixer..."
	@ $(CXX) -o $(BIN_DIR)/mesh_fixer $(FIXER_OBJS) $(LDFLAGS) $(CXX_LIBS)

implicit_uniform : $(SAMPLER_OBJS)
	@ echo "creating uniform implicit sampler..."
	@ $(CXX) -o $(BIN_DIR)/implicit_uniform $(SAMPLER_OBJS) $(LDFLAGS) $(CXX_LIBS)

single : $(SINGLE_OBJS)
	@ echo "creating single distribution plotter..."
	@ $(CXX) -o $(BIN_DIR)/single_distribution $(SINGLE_OBJS) $(LDFLAGS) $(CXX_LIBS)

aggregate : $(AGGREGATE_OBJS)
	@ echo "creating aggregate distribution plotter..."
	@ $(CXX) -o $(BIN_DIR)/aggregate_distribution $(AGGREGATE_OBJS) $(LDFLAGS) $(CXX_LIBS)

lib : $(LIB_OBJS)
	@ echo "Creating static library for evaluator ..."
	@ ar cr $(LIB_DIR)/libevaluator.a $(LIB_OBJS)
		@ ranlib $(LIB_DIR)/libevaluator.a

$(LIB_DIR)/%.o : %.cpp
	@ echo "Compiling $< ..."
	@ $(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(LIB_DIR)/*.o $(BIN_DIR)/run_evaluation $(BIN_DIR)/implicit_uniform $(BIN_DIR)/single_distribution $(BIN_DIR)/aggregate_distribution
