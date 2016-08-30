BABYDIR := txt/variables
LIBFILE := $(OBJDIR)/libStatObj.a

vpath %.cpp $(SRCDIR)
vpath %.cxx $(SRCDIR)
vpath %.hpp $(INCDIR)
vpath %.o $(OBJDIR)
vpath %.exe $(EXEDIR)
vpath %.d $(MAKEDIR)

CXX := $(shell root-config --cxx)
EXTRA_WARNINGS := -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Winit-self -Winvalid-pch -Wlong-long -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wpointer-arith -Wredundant-decls -Wstack-protector -Wswitch-default -Wswitch-enum -Wundef -Wunused -Wvariadic-macros -Wwrite-strings -Wabi -Wctor-dtor-privacy -Wnon-virtual-dtor -Wsign-promo -Wsign-compare #-Wunsafe-loop-optimizations -Wfloat-equal -Wsign-conversion -Wunreachable-code
CXXFLAGS := -isystem $(shell root-config --incdir) -Wall -Wextra -pedantic -Werror -Wshadow -Woverloaded-virtual -Wold-style-cast $(EXTRA_WARNINGS) $(shell root-config --cflags) -O2 -I $(INCDIR) -std=c++11
LD := $(shell root-config --ld)
LDFLAGS := $(shell root-config --ldflags)
LDLIBS := $(shell root-config --libs) -lMinuit -lRooStats -lRooFitCore -lRooFit -lTreePlayer 

GET_DEPS = $(CXX) $(CXXFLAGS) -MM -MP -MT "$(subst $(SRCDIR),$(OBJDIR),$(subst .cxx,.o,$(subst .cpp,.o,$<))) $@" -MF $@ $<
COMPILE = $(CXX) $(CXXFLAGS) -o $@ -c $<
LINK = $(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

BABY_FILES := $(wildcard $(BABYDIR)/*)
BABY_TYPES := $(notdir $(basename $(BABY_FILES)))
BABY_SRCS := $(addprefix $(SRCDIR)/core/baby_, $(addsuffix .cpp, $(BABY_TYPES)))
BABY_INCS := $(addprefix $(INCDIR)/core/baby_, $(addsuffix .hpp, $(BABY_TYPES)))
BABY_OBJS := $(addprefix $(OBJDIR)/core/baby_, $(addsuffix .o, $(BABY_TYPES)))
BABY_DEPS := $(addprefix $(MAKEDIR)/core/baby_, $(addsuffix .d, $(BABY_TYPES)))

FILTER_OUT = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))

HEADERS := $(call FILTER_OUT,.\#,$(shell find $(INCDIR) -name "*.hpp"))
OBJSRCS := $(call FILTER_OUT,.\#,$(shell find $(SRCDIR) -name "*.cpp"))
EXESRCS := $(call FILTER_OUT,.\#,$(shell find $(SRCDIR) -name "*.cxx"))
ALLSRCS := $(OBJSRCS) $(EXESRCS)

EXECUTABLES := $(subst $(SRCDIR),$(EXEDIR),$(subst .cxx,.exe,$(EXESRCS)))
OBJECTS := $(subst $(SRCDIR),$(OBJDIR),$(subst .cpp,.o,$(OBJSRCS))) $(OBJDIR)/core/baby.o $(BABY_OBJS)
DEPFILES := $(subst $(SRCDIR),$(MAKEDIR),$(subst .cpp,.d,$(subst .cxx,.d,$(ALLSRCS))))

PRINT_FUNC = echo -e "\e[34;1m$(1):\e[0m $($(1))"

all: delay

print_vars:
	$(call PRINT_FUNC,SRCDIR)
	$(call PRINT_FUNC,INCDIR)
	$(call PRINT_FUNC,MAKEDIR)
	$(call PRINT_FUNC,OBJDIR)
	$(call PRINT_FUNC,EXEDIR)

	$(call PRINT_FUNC,BABY_FILES)
	$(call PRINT_FUNC,BABY_TYPES)
	$(call PRINT_FUNC,BABY_SRCS)
	$(call PRINT_FUNC,BABY_INCS)
	$(call PRINT_FUNC,BABY_OBJS)
	$(call PRINT_FUNC,BABY_DEPS)

	$(call PRINT_FUNC,HEADERS)
	$(call PRINT_FUNC,OBJSRCS)
	$(call PRINT_FUNC,EXESRCS)
	$(call PRINT_FUNC,ALLSRCS)

	$(call PRINT_FUNC,EXECUTABLES)
	$(call PRINT_FUNC,OBJECTS)
	$(call PRINT_FUNC,DEPFILES)

include .subdirs.mk

.subdirs.mk: $(BABY_SRCS) $(BABY_INCS) $(HEADERS) $(ALLSRCS)
	./compile.py set_dirs --src_dir $(SRCDIR) --inc_dir $(INCDIR) --obj_dir $(OBJDIR) --make_dir $(MAKEDIR) --exe_dir $(EXEDIR)

$(OBJDIR)/core/generate_baby.o: $(SRCDIR)/core/generate_baby.cxx
	$(COMPILE)

$(EXEDIR)/core/generate_baby.exe: $(OBJDIR)/core/generate_baby.o
	$(LINK)

.SECONDARY: dummy_baby.all
.PRECIOUS: generate_baby.o

$(BABY_SRCS) $(BABY_INCS) $(SRCDIR)/core/baby.cpp $(INCDIR)/core/baby.hpp: dummy_baby.all
	: "Regenerated baby source code $@"

dummy_baby.all: $(EXEDIR)/core/generate_baby.exe $(BABY_FILES) $(BABYDIR)
	rm -f src/core/baby*.cpp inc/core/baby*.hpp bin/core/baby*.o bin/core/baby*.d
	./$< $(BABY_TYPES)

include $(DEPFILES) $(BABY_DEPS)

delay: $(EXECUTABLES)

$(LIBFILE): $(OBJECTS)

$(OBJDIR)/%.a:
	ar rcsv $@ $^

.PHONY: all delay test
.DELETE_ON_ERROR:
