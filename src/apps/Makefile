APPS	= leandvb leansdrscan
APPS	+= leansdrcat leantsgen leanchansim leandvbtx

all:	$(APPS)

clean:
	rm -f $(APPS)

DEPS	= ../leansdr/*.h

CXXFLAGS	= -O3 -I.. -Wall -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable

%:	%.cc $(DEPS)
	g++ $(CXXFLAGS) -DGUI $< -lX11 -o $@  ||  \
	g++ $(CXXFLAGS)       $<       -o $@

leandvb.embedded:	leandvb.cc $(DEPS)
	g++ $(CXXFLAGS) $< -static -o $@  ||  \
	g++ $(CXXFLAGS) $<         -o $@
