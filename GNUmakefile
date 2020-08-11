CXX      := g++
CXXFLAGS += -I. $(shell root-config --cflags) -g -I./Fade3D/include_fade3d -L./Fade3D/lib_ubuntu17.04_x86_64
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g -lfade3d -Wl,-rpath=./Fade3D/lib_ubuntu17.04_x86_64

PROGRAMS = CalibSCE MakeSmoothHists MakeSystVar

all:		clean $(PROGRAMS)

$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cpp -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
