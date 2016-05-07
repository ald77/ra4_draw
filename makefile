EXEDIR := run
OBJDIR := bin
SRCDIR := src
INCDIR := inc
MAKEDIR := bin
BABYDIR := txt/variables
LIBFILE := $(OBJDIR)/libStatObj.a

CXX := $(shell root-config --cxx)
EXTRA_WARNINGS := -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Winit-self -Winvalid-pch -Wlong-long -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wpointer-arith -Wredundant-decls -Wstack-protector -Wswitch-default -Wswitch-enum -Wundef -Wunused -Wvariadic-macros -Wwrite-strings -Wabi -Wctor-dtor-privacy -Wnon-virtual-dtor -Wsign-promo -Wsign-compare #-Wunsafe-loop-optimizations -Wfloat-equal -Wsign-conversion -Wunreachable-code
CXXFLAGS := -isystem $(shell root-config --incdir) -Wall -Wextra -pedantic -Werror -Wshadow -Woverloaded-virtual -Wold-style-cast $(EXTRA_WARNINGS) $(shell root-config --cflags) -O2 -I $(INCDIR) -std=c++11
LD := $(shell root-config --ld)
LDFLAGS := $(shell root-config --ldflags)
LDLIBS := $(shell root-config --libs) -lMinuit -lRooStats -lRooFitCore -lRooFit -lTreePlayer 

BABY_FILES := $(wildcard $(BABYDIR)/*)
BABY_TYPES := $(notdir $(basename $(BABY_FILES)))
BABY_SRCS := $(addprefix $(SRCDIR)/baby_, $(addsuffix .cpp, $(BABY_TYPES)))
BABY_INCS := $(addprefix $(INCDIR)/baby_, $(addsuffix .hpp, $(BABY_TYPES)))
BABY_OBJS := $(addprefix $(OBJDIR)/baby_, $(addsuffix .o, $(BABY_TYPES)))
BABY_DEPS := $(addprefix $(MAKEDIR)/baby_, $(addsuffix .d, $(BABY_TYPES)))

EXECUTABLES := $(addprefix $(EXEDIR)/, $(addsuffix .exe, $(notdir $(basename $(wildcard $(SRCDIR)/*.cxx))))) 
OBJECTS := $(addprefix $(OBJDIR)/, $(addsuffix .o, $(notdir $(basename $(wildcard $(SRCDIR)/*.cpp))))) $(OBJDIR)/baby.o $(BABY_OBJS)

FIND_DEPS = $(CXX) $(CXXFLAGS) -MM -MG -MF $@ $<
EXPAND_DEPS = perl -pi -e 's|$*.o|$(OBJDIR)/$*.o $(MAKEDIR)/$*.d|g' $@
GET_DEPS = $(FIND_DEPS) && $(EXPAND_DEPS)
COMPILE = $(CXX) $(CXXFLAGS) -o $@ -c $<
LINK = $(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

vpath %.cpp $(SRCDIR)
vpath %.cxx $(SRCDIR)
vpath %.hpp $(INCDIR)
vpath %.o $(OBJDIR)
vpath %.exe $(EXEDIR)
vpath %.d $(MAKEDIR)

all: $(EXECUTABLES)

-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cpp)))))
-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cxx)))))
-include $(MAKEDIR)/baby.d $(BABY_DEPS)

$(LIBFILE): $(OBJECTS)

$(MAKEDIR)/%.d: $(SRCDIR)/%.cpp
	$(GET_DEPS)

$(MAKEDIR)/%.d: $(SRCDIR)/%.cxx
	$(GET_DEPS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILE)

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	$(COMPILE)

$(OBJDIR)/%.a:
	ar rcsv $@ $^

$(EXEDIR)/generate_baby.exe: $(OBJDIR)/generate_baby.o
	$(LINK)

$(EXEDIR)/%.exe: $(OBJDIR)/%.o $(LIBFILE)
	$(LINK)

.SECONDARY: dummy_baby.all
.PRECIOUS: generate_baby.o

$(BABY_SRCS) $(BABY_INCS) $(SRCDIR)/baby.cpp $(INCDIR)/baby.hpp: dummy_baby.all
	echo "Regenerated baby source code $@"

dummy_baby.all: $(EXEDIR)/generate_baby.exe $(BABY_FILES) $(BABYDIR)
	rm -f src/baby*.cpp inc/baby*.hpp bin/baby*.o bin/baby*.d
	./$< $(BABY_TYPES)

.DELETE_ON_ERROR:
