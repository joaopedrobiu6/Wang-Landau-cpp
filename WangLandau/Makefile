#######################################################################
#### Definitions
#######################################################################
BINDIR := bin
LIBDIR := lib
MAINDIR := main
SRCDIR := src
RESUL := results

CXX := g++
CXXFLAGS := -O3 -std=c++23 -I $(SRCDIR) -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include

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
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(LIBDIR) -l CompPhys -L/opt/homebrew/opt/libomp/lib -lomp

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp | $(LIBDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(BINDIR) $(LIBDIR)

.PHONY: all clean
