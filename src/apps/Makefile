default:	native

VERSION	:= leansdr-$(shell git describe)

# Portable flags for development.

CXXFLAGS	= -O3 -I.. -DVERSION=\"$(VERSION)\" \
		  -Wall \
		  -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable

LDFLAGS		=

# Base applications have no external dependencies.
APPS	=
APPS	+= leandvbtx leandvb
APPS	+= leansdrscan leansdrserv
APPS	+= leantsgen leansdrcat leanchansim


# leanmlmrx requires fftw3 and pthread.
ifneq "$(wildcard /usr/include/fftw3.h)" ""
ifneq "$(wildcard /usr/include/pthread.h)" ""
APPS		+= leanmlmrx
CXXFLAGS	+=
LDFLAGS		+= -lfftw3f -lpthread
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
CXXFLAGS	+= -Ilibiio
LDFLAGS		+= -Llibiio -L. -liio -lpthread -lxml2 -lz -llzma
else
$(info libiio not found. Will not build leaniio{rx,tx}.)
endif

# Enable X11 GUI if easily found.
ifneq "$(wildcard /usr/include/X11/Xlib.h)" ""
CXXFLAGS_GUI	= -DGUI
LDFLAGS_GUI	= -lX11
else
$(info X11 not found. Will not compile GUI.)
endif

DEPS	= ../leansdr/*.h

%:	%.cc $(DEPS)
	g++ $(CXXFLAGS) $(CXXFLAGS_GUI) $< $(LDFLAGS) $(LDFLAGS_GUI) -o $@

native::	$(APPS)
	@echo "Now 'make embedded' to also build static binaries."

clean::
	rm -f $(APPS)

##### CPU-specific flags for embedded platforms.

ARCH	:= $(shell uname -m)

ifneq "$(findstring arm,$(ARCH))" ""
CXXFLAGS_EMBEDDED	= $(CXXFLAGS) -Ofast -mfpu=neon -funsafe-math-optimizations -fsingle-precision-constant
LDFLAGS_EMBEDDED	= $(LDFLAGS) -static
else
CXXFLAGS_EMBEDDED	= $(CXXFLAGS)
LDFLAGS_EMBEDDED	= $(LDFLAGS) -static
endif

embedded:	$(APPS:%=%.$(ARCH))

%.$(ARCH):	%.cc $(DEPS)
	g++ $(CXXFLAGS_EMBEDDED) $< $(LDFLAGS_EMBEDDED) -o $@
