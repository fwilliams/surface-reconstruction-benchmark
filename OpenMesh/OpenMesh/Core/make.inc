SRC_DIRS:=$(SRC_DIRS):$(PARENT_DIR)/Core

SRCDIR=$(PARENT_DIR)/Core

ALL_OBJS := $(ALL_OBJS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.cc))
ALL_DEPS := $(ALL_DEPS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.d,$(wildcard $(SRCDIR)/*.cc))

PARENT_DIR=OpenMesh/Core
include $(PARENT_DIR)/Geometry/make.inc

PARENT_DIR=OpenMesh/Core
include $(PARENT_DIR)/IO/make.inc

PARENT_DIR=OpenMesh/Core
include $(PARENT_DIR)/Mesh/make.inc

PARENT_DIR=OpenMesh/Core
include $(PARENT_DIR)/System/make.inc

PARENT_DIR=OpenMesh/Core
include $(PARENT_DIR)/Templates/make.inc

PARENT_DIR=OpenMesh/Core
include $(PARENT_DIR)/Utils/make.inc
