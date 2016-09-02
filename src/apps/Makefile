APPS	= leandvb leansdrscan leansdrcat leandvb_gui 

all:	$(APPS)

DEPS	= ../leansdr/*.h

CXXFLAGS	= -O3 -I.. -Wall -Wno-sign-compare -Wno-array-bounds -Wno-unused-variable

%:	%.cc $(DEPS)
	g++ $(CXXFLAGS) $< -o $@

%_gui:	%.cc $(DEPS)
	g++ $(CXXFLAGS) -DGUI $< -lX11 -o $@
