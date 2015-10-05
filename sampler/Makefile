MODELING_LIB = ../modeling/lib

TEEM_FLAG=0

CXX=g++

OPENMESH_INC = ../OpenMesh
ANN_INC = ../ann/include

OPENMESH_LIB = ../OpenMesh/lib
ANN_LIB = ../ann/lib

CXXFLAGS = -g -I$(OPENMESH_INC) -I$(ANN_INC) -I/usr/include/libpng14
LDFLAGS = -lblas -llapack -L$(MODELING_LIB) -lmodeling -L$(OPENMESH_LIB) -lOpenMeshCore -L$(ANN_LIB) -lANN -lpng14

CXXFLAGS += -DHAS_TEEM=0

ifeq ($(TERM_PROGRAM), Apple_Terminal)
	CXXFLAGS += -DHAS_MALLOC=0
else
	CXXFLAGS += -DHAS_MALLOC=1
endif

LIB_DIR=lib
BIN_DIR=../bin

SAMPLER_OBJS = $(LIB_DIR)/ImplicitNewton.o $(LIB_DIR)/UniformSampler.o \
					$(LIB_DIR)/RangeImage.o $(LIB_DIR)/OrientedRangeImage.o \
					$(LIB_DIR)/RayCaster.o $(LIB_DIR)/UniformOccupancyGrid.o $(LIB_DIR)/RangeScanner.o \
					$(LIB_DIR)/GaussianFit.o $(LIB_DIR)/PNGImage.o $(LIB_DIR)/SphericalImager.o

PC_OBJS = $(LIB_DIR)/PointCloud.o $(LIB_DIR)/OrientedPointCloud.o $(LIB_DIR)/XYZReader.o $(LIB_DIR)/NPTSReader.o \
			 $(LIB_DIR)/SparseRangeScan.o

RUN_UNIFORM_OBJ = $(SAMPLER_OBJS) $(PC_OBJS) $(LIB_DIR)/run_uniform.o
RUN_PC_GENERATOR_OBJ = $(SAMPLER_OBJS) $(PC_OBJS) $(LIB_DIR)/pc_generator.o
POINT_ORIENTER_OBJ = $(SAMPLER_OBJS) $(PC_OBJS) $(LIB_DIR)/point_orienter.o
PROJECTION_OBJ = $(SAMPLER_OBJS) $(PC_OBJS) $(LIB_DIR)/projection_cleanup.o

LIB_OBJS = $(SAMPLER_OBJS) $(PC_OBJS)

default : all
all : lib uniform generator orienter point_projection

uniform : $(RUN_UNIFORM_OBJ)
	@ echo "creating uniform"
	@ $(CXX) -o $(BIN_DIR)/uniform_sampler $(RUN_UNIFORM_OBJ) $(LDFLAGS)

generator : $(RUN_PC_GENERATOR_OBJ)
	@ echo "creating pc configuration generator"
	@ $(CXX) -o $(BIN_DIR)/pc_generator $(RUN_PC_GENERATOR_OBJ) $(LDFLAGS)

orienter : $(POINT_ORIENTER_OBJ)
	@ echo "creating point orienter"
	@ $(CXX) -o $(BIN_DIR)/point_orienter $(POINT_ORIENTER_OBJ) $(LDFLAGS)

point_projection : $(PROJECTION_OBJ)
	@ echo "creating point projection"
	@ $(CXX) -o $(BIN_DIR)/point_projection $(PROJECTION_OBJ) $(LDFLAGS)

lib : $(LIB_OBJS)
	@ echo "Creating static library for sampler ..."
	@ ar cr $(LIB_DIR)/libsampler.a $(LIB_OBJS)
		@ ranlib $(LIB_DIR)/libsampler.a

$(LIB_DIR)/%.o : %.cpp
	@ echo "Compiling $< ..."
	@ $(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	rm -f $(LIB_DIR)/*.o $(LIB_DIR)/libsampler.a $(BIN_DIR)/uniform_sampler $(BIN_DIR)/pc_generator $(BIN_DIR)/point_orienter $(BIN_DIR)/point_projection
