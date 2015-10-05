SRC_DIRS:=$(SRC_DIRS):$(PARENT_DIR)/OpenMesh

SRCDIR=$(PARENT_DIR)/OpenMesh

ALL_OBJS := $(ALL_OBJS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.cc))
ALL_DEPS := $(ALL_DEPS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.d,$(wildcard $(SRCDIR)/*.cc))

PARENT_DIR=OpenMesh
include $(PARENT_DIR)/Core/make.inc

PARENT_DIR=OpenMesh
include $(PARENT_DIR)/Tools/make.inc
