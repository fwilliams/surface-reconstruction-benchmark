# Benedict Brown
# Princeton University
#
# Makefile
# Makefile for tps alignment
# Based on the trimesh2 Makefile

MAKERULESDIR = ..

include $(MAKERULESDIR)/Makerules

BINDIR = ../../../bin

COMMON_OBJS = $(addprefix $(OBJDIR)/, KDtree.o ICP.o EdgeMesh.o)
DEPLIBS = ../../trimesh2/lib.$(TRIMESH_UNAME)/libtrimesh.a
#INCLUDES += -I../include -I../../glui-2.35/src/include/GL -I../../trimesh2/include -I../../tnt
INCLUDES += -I../include -I../../trimesh2/include/GL -I../../trimesh2/include -I../../tnt
LIBS := $(COMMON_OBJS) $(LIBS)

MAKEFILE_DEPS := Makefile ../Makerules ../Makedefs.$(UNAME)

EXES     = autorient correspond preprocess make_session distance pre_to_ply ply2vn
ALL_EXES = $(EXES) interactive global_reg

default: $(BINDIR) $(addprefix $(BINDIR)/, $(ALL_EXES))

$(BINDIR):
	@ -mkdir $@

clean:
	rm -rf $(OBJDIR)/* \
	rm $(BINDIR)/autorient $(BINDIR)/correspond $(BINDIR)/preprocess $(BINDIR)/make_session $(BINDIR)/distance $(BINDIR)/pre_to_ply $(BINDIR)/ply2vn $(BINDIR)/interactive $(BINDIR)/global_reg

# clobber: clean
# 	rm -f *~

$(BINDIR)/global_reg : $(OBJDIR)/tps.o $(OBJDIR)/global_reg_util.o
	$(LINK) $(OBJDIR)/tps.o $(OBJDIR)/global_reg_util.o

$(BINDIR)/interactive : $(OBJDIR)/tps.o
	$(LINK) $(OBJDIR)/tps.o $(XLIBS) $(GLUI_LIBS)

$(addprefix $(BINDIR)/, $(EXES)) :
	$(LINK)

#%.o:%.cc Makefile
#	$(CXX) $(CFLAGS) $(INCLUDES) -c -o $@ $<

#%.o:%.c Makefile
#	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# dependencies
$(addprefix $(BINDIR)/, $(ALL_EXES)): $(BINDIR)/%: $(BINDIR) $(COMMON_OBJS) $(DEPLIBS) $(MAKEFILE_DEPS) $(OBJDIR)/%.o
