SRC_DIRS:=$(SRC_DIRS):$(PARENT_DIR)/Tools

SRCDIR=$(PARENT_DIR)/Tools

ALL_OBJS := $(ALL_OBJS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.cc))
ALL_DEPS := $(ALL_DEPS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.d,$(wildcard $(SRCDIR)/*.cc))

PARENT_DIR=OpenMesh/Tools
include $(PARENT_DIR)/Utils/make.inc
