CC = gcc
CXX = g++
AR = ar

COPTS += -O3 -ffast-math -funroll-loops -fomit-frame-pointer
COPTS += -Wall -W -Wno-unused -Wno-sign-compare -Wno-unknown-pragmas
CXXOPTS = $(COPTS)

CXXFLAGS = -I./

OBJDIR = lib
PARENT_DIR =.
include OpenMesh/make.inc
VPATH=$(SRC_DIRS)

$(OBJDIR)/%.o : %.cc
	@ echo "Compiling $< ..."
	@ $(CXX) -c $(CXXFLAGS) $< -o $@

openmesh : $(ALL_OBJS)
	@ echo "Creating static library for OpenMeshCore..."
	@ ar cr $(OBJDIR)/libOpenMeshCore.a $(ALL_OBJS)
		@ ranlib $(OBJDIR)/libOpenMeshCore.a

-include $(ALL_DEPS)

$(OBJDIR)/%.d: %.cc
	@echo "Building dependencies for $< ..."
	@set -e; rm -f $@; \
		$(CC) -MM $(CXXFLAGS) $< > $(OBJDIR)/$(*F).t1; \
		sed 's,\(.*\)\.o[ :]*,$(OBJDIR)/$(*F).o $(OBJDIR)/$(@F) : ,g' < $(OBJDIR)/$(*F).t1 > $(OBJDIR)/$(*F).t2; \
		sed '1n; /.*\.cc/d' < $(OBJDIR)/$(*F).t2 > $(OBJDIR)/$(@F); \
		rm -f $(OBJDIR)/$(*F).t1 $(OBJDIR)/$(*F).t2

.PHONY : clean
clean:
	@ echo "cleaning object files and library ..."
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d $(OBJDIR)/libOpenMeshCore.a
