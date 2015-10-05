SRC_DIRS:=$(SRC_DIRS):$(PARENT_DIR)/Utils

SRCDIR=$(PARENT_DIR)/Utils

ALL_OBJS := $(ALL_OBJS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.cc))
ALL_DEPS := $(ALL_DEPS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.d,$(wildcard $(SRCDIR)/*.cc))
