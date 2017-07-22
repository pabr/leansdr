APPS	= leandvb leansdrscan
APPS	+= leansdrcat leantsgen leanchansim leandvbtx

all:	$(APPS)

clean:
	rm -f $(APPS)

DEPS	= ../leansdr/*.h

CXXFLAGS = -O3 -I.. -Wall -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable -Wno-tautological-compare
CXXFLAGS += -L/opt/X11/lib -lX11 -lpthread

%:	%.cc $(DEPS)
	g++ $(CXXFLAGS) -DGUI $< -lX11 -o $@  ||  \
	g++ $(CXXFLAGS)       $<      -o $@

EMBED_FLAGS=	-I.. -Wall -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable  \
		-Ofast -lm -mfpu=neon -funsafe-math-optimizations -fsingle-precision-constant -lm

leandvb.embedded:	leandvb.cc $(DEPS)
	clang $(CXXFLAGS) $(EMBED_FLAGS) $< -static -o $@  ||  \
	clang $(CXXFLAGS) $(EMBED_FLAGS) $<         -o $@
