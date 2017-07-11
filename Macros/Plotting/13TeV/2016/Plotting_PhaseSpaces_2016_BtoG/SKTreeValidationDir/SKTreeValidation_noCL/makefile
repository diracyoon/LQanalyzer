all: MakeDataMCCompPlots

OBJDIR = obj
DEPDIR  = $(OBJDIR)/dep
EXECPP = MakeDataMCCompPlots.cpp
LIBRARY=Plotting

ARCH_LOC_1 := $(wildcard $(shell root-config --prefix)/test/Makefile.arch)
ARCH_LOC_2 := $(wildcard $(shell root-config --prefix)/share/root/test/Makefile.arch)
ARCH_LOC_3 := $(wildcard $(shell root-config --prefix)/share/doc/root/test/Makefile.arch)
ARCH_LOC_4 := $(wildcard $(shell root-config --prefix)/etc/Makefile.arch)
ARCH_LOC_5 := $(wildcard $(shell root-config --prefix)/etc/root/Makefile.arch)
ARCH_LOC_5 := $(wildcard $(shell root-config --prefix)/etc/root/Makefile.arch)
ARCH_LOC_6 := $(wildcard $(shell root-config --etcdir)/Makefile.arch)
ifneq ($(strip $(ARCH_LOC_1)),)
  $(info Using $(ARCH_LOC_1))
  include $(ARCH_LOC_1)
else
  ifneq ($(strip $(ARCH_LOC_2)),)
    $(info Using $(ARCH_LOC_2))
    include $(ARCH_LOC_2)
  else
    ifneq ($(strip $(ARCH_LOC_3)),)
      $(info Using $(ARCH_LOC_3))
      include $(ARCH_LOC_3)
    else
                ifneq ($(strip $(ARCH_LOC_4)),)
        $(info Using $(ARCH_LOC_4))
        include $(ARCH_LOC_4)
      else
        ifneq ($(strip $(ARCH_LOC_5)),)
          $(info Using $(ARCH_LOC_5))
          include $(ARCH_LOC_5)
        else
         ifneq ($(strip $(ARCH_LOC_6)),)
          $(info Using $(ARCH_LOC_6))
          include $(ARCH_LOC_6)
         else
         $(error Could not find Makefile.arch!)
        endif
        endif
        endif
    endif
  endif
endif

# compilation options
VPATH    += $(OBJDIR)
INCLUDES += -I./ -I$(MACRO_DIR)
CXXFLAGS += -Wall -Wno-overloaded-virtual -Wno-unused
SRCSUF = cpp
ROOTLIBS = `root-config --libs`

CPPLIST = $(filter-out $(EXECPP), $(wildcard *.$(SRCSUF)) )
OLIST = $(patsubst %.$(SRCSUF),%.o,$(notdir $(CPPLIST)))

LIBFILE   = ./lib$(LIBRARY).a
SHLIBFILE = ./lib$(LIBRARY).$(DllSuf)

# Implicit rule to compile all sources
%.o : %.$(SRCSUF)
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CXX) $(CXXFLAGS) -O2 -c $< -o $(OBJDIR)/$(notdir $@) $(INCLUDES)


##############################
# The dependencies section   
# - the purpose of the .d files is to keep track of the
#   header file dependence
# - this can be achieved using the makedepend command 
##############################
# .d tries to pre-process .cc
-include $(foreach var,$(notdir $(CPPLIST:.$(SRCSUF)=.d)),$(DEPDIR)/$(var))

$(DEPDIR)/%.d: %.$(SRCSUF)
	@mkdir -p $(DEPDIR)
	if test -f $< ; then \
		echo "Making $(@F)"; \
		$(SHELL) -ec '$(CPP) -MM $(CXXFLAGS) $(INCLUDES) $< | sed '\''/Cstd\/rw/d'\'' > $@'; \
	fi

# Rule to combine objects into a unix shared library
$(SHLIBFILE): $(OLIST)
	@echo "Making shared library: $(SHLIBFILE)"
	@rm -f $(SHLIBFILE)
ifneq (,$(findstring macosx,$(ARCH)))
	@$(LD) $(LDFLAGS) -dynamiclib -single_module -undefined dynamic_lookup -O2 $(addprefix $(OBJDIR)/,$(OLIST)) -o $(SHLIBFILE)
else
	@$(LD) $(LDFLAGS) $(SOFLAGS) -O2 $(addprefix $(OBJDIR)/,$(OLIST)) -o $(SHLIBFILE)
endif

# Rule to make the executable
MakeDataMCCompPlots: $(SHLIBFILE) MakeDataMCCompPlots.cpp
	@echo "Make executable: $@"
	@$(CXX) $(CXXFLAGS) -lASImage  -o $@  $(SHLIBFILE) MakeDataMCCompPlots.cpp $(ROOTLIBS) $(INCLUDES)


clean:
	@rm -f $(DICTFILE) $(DICTHEAD)
	@rm -f $(OBJDIR)/*.o
	@rm -f $(SHLIBFILE)
	@rm -f $(SFRAME_LIB_PATH)/lib$(LIBRARY).so

distclean:
	@rm -rf $(OBJDIR)
	@rm -rf $(DEPDIR)
	@rm -f $(DICTFILE) $(DICTHEAD)
	@rm -f $(SHLIBFILE)