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

EMBED_FLAGS=	-I.. -Wall -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable  \
		-Ofast -mfpu=neon -funsafe-math-optimizations -fsingle-precision-constant

leandvb.embedded:	leandvb.cc $(DEPS)
	g++ $(CXXFLAGS) $(EMBED_FLAGS) $< -static -o $@  ||  \
	g++ $(CXXFLAGS) $(EMBED_FLAGS) $<         -o $@
