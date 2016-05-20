CPPFLAGS = -I include -Wall -std=c++11
CPPFLAGS += -O2
LDLIBS = -lm -lrt


FILES   = background_subtraction


INCLUDES = 
LIBRARY += $(LIBDIR)/descriptor -ldescriptor
LIBRARY += $(LIBDIR)/ImageContent -limage
LIBRARY += $(LIBDIR)/segment -lsegment
LIBRARY += $(LIBDIR)/gauss_iir -lgauss_iir
LIBRARY += $(LIBDIR)/util -ltools
LIBRARY += $(LIBDIR)/matrix -lmatrix
LIBRARY += $(LIBDIR)/util -ltools
LIBRARY += -L/opt/local/lib	   -lpng  -ljpeg -lz -lm 

include ../makefile
