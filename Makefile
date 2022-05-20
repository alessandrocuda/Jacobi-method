UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXX = g++
endif
ifeq ($(UNAME_S),Darwin)
	FF_LOCAL = -I $$FF_ROOT
    CXX = g++-11
endif

# define C++ compiler to use

# define any compile-time flags
CFLAGS = -std=c++11 -O3 -funsafe-math-optimizations #-Wall -g

LDFLAGS = -lpthread #-lboost_system -lcrypto -lssl -lcpprest -lpthread

OBJDIR = build
SRCDIR = src
TESTDIR = test

SRC := $(shell find $(SRCDIR) -name "*.cpp")
SRC_OBJS = $(SRC:%.cpp=$(OBJDIR)/%.o)

JACOBI_TEST_SRC	:= $(TESTDIR)/test.cpp
JACOBI_TEST_OBJS = $(JACOBI_TEST_SRC:%.cpp=$(OBJDIR)/%.o)
JACOBI_TEST = bin/test

PAR_ANALYSIS_SRC	:= $(TESTDIR)/parallel_analysis.cpp
PAR_ANALYSIS_OBJS = $(PAR_ANALYSIS_SRC:%.cpp=$(OBJDIR)/%.o)
PAR_ANALYSIS = bin/parallel_analysis

all:    $(JACOBI_TEST) $(PAR_ANALYSIS)
		@echo parallel_analysis, test have been compiled.

$(JACOBI_TEST): $(SRC_OBJS) $(JACOBI_TEST_OBJS)
	@echo "== LINKING EXECUTABLE $@"
	@mkdir -p $(@D)
	@$(CXX) $^ $(LDFLAGS) $(FF_LOCAL) -o $(JACOBI_TEST)

$(PAR_ANALYSIS): $(SRC_OBJS) $(PAR_ANALYSIS_OBJS)
	@echo "== LINKING EXECUTABLE $@"
	@mkdir -p $(@D)
	@$(CXX) $^ $(LDFLAGS) $(FF_LOCAL) -o $(PAR_ANALYSIS)
	
$(OBJDIR)/%.o: %.cpp
	@echo "COMPILING SOURCE $< INTO OBJECT $@"
	@mkdir -p '$(@D)'
	@$(CXX) -c $(CXXFLAGS) $(FF_LOCAL) $(LDFLAGS) $< -o $@


clean:
	find . -name *.o -delete -print
	rm -f -v $(JACOBI_TEST) $(PAR_ANALYSIS)