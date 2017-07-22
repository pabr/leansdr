APPS	= leandvb leansdrscan
APPS	+= leansdrcat leantsgen leanchansim leandvbtx
APPS	+= leandvb.embedded 

all:	$(APPS)

clean:
	rm -f $(APPS)

DEPS	= ../leansdr/*.h

CXXFLAGS	= -O3 -I.. -Wall -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable
CXXFLAGS	+= -L/opt/X11/lib -lX11 -lpthread

%:	%.cc $(DEPS)
	g++ $(CXXFLAGS) -DGUI $< -lX11 -o $@  ||  \
	g++ $(CXXFLAGS)       $<       -o $@

EMBED_FLAGS=	-I.. -Wall -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable  \
		-Ofast -mfpu=neon -funsafe-math-optimizations -fsingle-precision-constant

leandvb.embedded:	leandvb.cc $(DEPS)
	g++ $(CXXFLAGS) $(EMBED_FLAGS) $< -static -o $@  ||  \
	g++ $(CXXFLAGS) $(EMBED_FLAGS) $<         -o $@
