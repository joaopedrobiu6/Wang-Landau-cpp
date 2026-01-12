#######################################################################
#### Definitions
#######################################################################
BINDIR := bin
LIBDIR := lib
MAINDIR := main
SRCDIR := src
RESUL := results

CXX := g++

# Detect libomp path
OMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
ifeq ($(OMP_PREFIX),)
  $(error libomp not found. Please install it with 'brew install libomp')
endif

OMP_INCLUDE := -I$(OMP_PREFIX)/include
OMP_LIB := -L$(OMP_PREFIX)/lib -lomp

CXXFLAGS := -O3 -std=c++23 -I $(SRCDIR) -Xpreprocessor -fopenmp $(OMP_INCLUDE)

#######################################################################
#### Files
#######################################################################
SRCS_CPP := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(patsubst %.cpp, $(LIBDIR)/%.o, $(notdir $(SRCS_CPP)) )
MAIN := $(patsubst %.cpp, $(BINDIR)/%.exe, $(notdir $(wildcard $(MAINDIR)/*.cpp)) )

#######################################################################
#### Rules
#######################################################################
all: $(MAIN)

$(BINDIR) $(LIBDIR) $(RESUL):
	mkdir -p $@

$(LIBDIR)/libCompPhys.a: $(OBJS) | $(LIBDIR)
	ar ruv $@ $^
	ranlib $@
	rm -f $^

$(BINDIR)/%.exe: $(MAINDIR)/%.cpp $(LIBDIR)/libCompPhys.a | $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(LIBDIR) -l CompPhys $(OMP_LIB)

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp | $(LIBDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(BINDIR) $(LIBDIR)

.PHONY: all clean
