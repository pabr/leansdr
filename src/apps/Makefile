APPS	= leandvb leansdrscan leansdrcat

all:	$(APPS)

clean:
	rm -f $(APPS)

DEPS	= ../leansdr/*.h

CXXFLAGS	= -O3 -I.. -Wall -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable

%:	%.cc $(DEPS)
	g++ $(CXXFLAGS) -DGUI $< -lX11 -o $@  ||  \
	g++ $(CXXFLAGS)       $<       -o $@
