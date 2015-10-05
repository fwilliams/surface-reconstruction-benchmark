SRC_DIRS:=$(SRC_DIRS):$(PARENT_DIR)/reader

SRCDIR=$(PARENT_DIR)/reader

ALL_OBJS := $(ALL_OBJS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.cc))
ALL_DEPS := $(ALL_DEPS) $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.d,$(wildcard $(SRCDIR)/*.cc))
