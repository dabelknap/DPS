DIR=$(shell pwd)

LIBPATH = -L$(DIR)/pythia8170/lib/archive
HDRPATH = -I$(DIR)/pythia8170/include

gen_events : gen_events.cc
	g++ -o $@ $< $(LIBPATH) $(HDRPATH) -llhapdfdummy -lpythia8 `root-config --cflags --glibs`
