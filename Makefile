CXX := g++
CXXFLAGS := -DRESTRICT=__restrict__ -Drestrict=__restrict__ -export-dynamic -fPIC -std=c++11 -O3 -DNDEBUG
LDFLAGS := -lnlopt
INCLUDES := -Iinclude -Ipsi4-missing-headers/ -I/home/rmcgibbo/miniconda/envs/2.7/bin/../include/ -I/home/rmcgibbo/miniconda/envs/2.7/bin/../include/python2.7/

CXXDEFS := -DHAVE_DKH -DHAVE_MM_MALLOC_H -DHAVE_MKL_LAPACK -DHAVE_MKL_BLAS -DHAS_CXX11_VARIADIC_TEMPLATES -DHAS_CXX11_STATIC_ASSERT -DHAS_CXX11_SIZEOF_MEMBER -DHAS_CXX11_RVALUE_REFERENCES -DHAS_CXX11_NULLPTR -DHAS_CXX11_LONG_LONG -DHAS_CXX11_LAMBDA -DHAS_CXX11_INITIALIZER_LIST -DHAS_CXX11_DECLTYPE -DHAS_CXX11_CSTDINT_H -DHAS_CXX11_CONSTEXPR -DHAS_CXX11_AUTO_RET_TYPE -DHAS_CXX11_AUTO -DHAS_CXX11_FUNC -DHAS_CXX11 -DSYS_LINUX -DUSE_FCMANGLE_H -D_GLIBCXX_USE_CXX11_ABI=0

# The name of your plugin. Taken from the directory name.
NAME := $(shell basename `pwd`)

# Used to determine linking flags.
UNAME = $(shell uname)

PSITARGET = $(NAME).so

# Start the compilation rules
default:: $(PSITARGET)

# Add the flags needed for shared library creation
ifeq ($(UNAME), Darwin)
    CXXFLAGS += -fno-common
endif

CXXSRC := $(wildcard src/*.cc)
BINOBJ := $(addprefix obj/, $(notdir $(CXXSRC:%.cc=%.o)))

src/gitversion.hpp: .git/HEAD .git/index
	echo "const char *GIT_VERSION = \"$(shell git rev-parse --short HEAD)\";" > $@
obj/dotsphere.o: src/dotsphere.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
obj/potential.o: src/potential.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
obj/resp2.o: src/resp2.cc src/gitversion.hpp
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
obj/vdwsurface.o: src/vdwsurface.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
obj/respfit.o: src/respfit.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(PSITARGET): $(BINOBJ)
	$(CXX) -shared -o $@ $^ $(CXXDEFS) $(LDFLAGS)

clean: psiclean
	rm -f $(BINOBJ) $(PSITARGET) *.d *.pyc *.test output.dat psi.timer.dat

psiclean:
	rm -f `cat psi.*.clean | xargs`
	rm -f psi.*.clean  psi.*.clean.out timer.dat psi.timer.dat

print-%:
	@echo '$*=$($*)'
