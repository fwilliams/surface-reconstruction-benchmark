SRC_DIRS:=$(SRC_DIRS):$(PARENT_DIR)/IO

SRCDIR=$(PARENT_DIR)/IO

ALL_OBJS := $(ALL_OBJS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.cc))
ALL_DEPS := $(ALL_DEPS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.d,$(wildcard $(SRCDIR)/*.cc))

PARENT_DIR=OpenMesh/Core/IO
include $(PARENT_DIR)/exporter/make.inc

PARENT_DIR=OpenMesh/Core/IO
include $(PARENT_DIR)/importer/make.inc

PARENT_DIR=OpenMesh/Core/IO
include $(PARENT_DIR)/reader/make.inc

PARENT_DIR=OpenMesh/Core/IO
include $(PARENT_DIR)/writer/make.inc
