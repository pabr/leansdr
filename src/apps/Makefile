default:	generic

help:
	@echo "make generic      - With generic flags."
	@echo "make embedded     - Static, with platform-specific flags"
	@echo "make OUTDIR=...   - Change output directory"

OUTDIR	= .

VERSION	:= leansdr-$(shell git describe)

# Generic flags, for development.

CXXFLAGS	= -O3 -I.. -DVERSION=\"$(VERSION)\" \
		  -Wall \
		  -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable

LDFLAGS		=

DEPS	= ../leansdr/*.h

# Uncomment to compile features with potential intellectual property issues.
#CXXFLAGS	+= -DLEANSDR_EXTENSIONS

# Base applications have no external dependencies.
APPS	=
APPS	+= leandvbtx leandvb
APPS	+= leansdrscan leansdrserv
APPS	+= leantsgen leansdrcat leanchansim

ifneq "$(wildcard /usr/include/X11/Xlib.h)" ""
$(OUTDIR)/leandvb:	CXXFLAGS_APP=-DGUI
$(OUTDIR)/leandvb:	LDFLAGS_APP=-lX11
else
$(info libX11 not found. Will not support --gui.)
endif

# leanmlmrx requires fftw3 and pthread.
ifneq "$(wildcard /usr/include/fftw3.h)" ""
ifneq "$(wildcard /usr/include/pthread.h)" ""
APPS	+= leanmlmrx
$(OUTDIR)/leanmlmrx:	LDFLAGS_APP=-lfftw3f -lpthread
else
$(info libpthread not found. Will not build leanmlmrx.)
endif
else
$(info libfftw3f not found. Will not build leanmlmrx.)
endif

# leaniio{rx,tx} require libiio.
# Expected in ./libiio/, including libiio.a.
# (Assume libiio dependencies are satisfied.)
ifneq "$(wildcard libiio/iio.h)" ""
APPS		+= leaniiorx leaniiotx
$(OUTDIR)/leaniiorx $(OUTDIR)/leaniiotx:	CXXFLAGS_APP=-Ilibiio
$(OUTDIR)/leaniiorx $(OUTDIR)/leaniiotx:	LDFLAGS_APP=-Llibiio -Llibiio-$(ARCH) -L. -liio -lpthread -lxml2 -lz -llzma
else
$(info libiio not found. Will not build leaniio{rx,tx}.)
endif

$(OUTDIR)/%:	%.cc $(DEPS) $(DEPS_APP)
	g++ $(CXXFLAGS) $(CXXFLAGS_APP) $< $(LDFLAGS) $(LDFLAGS_APP) -o $@

generic:	$(APPS:%=$(OUTDIR)/%)

embedded:	$(APPS:%=$(OUTDIR)/%)
embedded:	LDFLAGS+=-static

clean:
	rm -f $(APPS:%=$(OUTDIR)/%)
