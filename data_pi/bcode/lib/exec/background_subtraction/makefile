CPPFLAGS = -I include -Wall -std=c++11
CPPFLAGS += -O2
CPPFLAGS += -march=armv7-a
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

#TBB libraries
TBB_DIR = /home/pi/tbb43_20150611oss
#CXXFLAGS += -DTBB_USE_GCC_BUILTINS=1 -D__TBB_64BIT_ATOMICS=0

TBB_INC_DIR ?= $(TBB_DIR)/include
TBB_LIB_DIR ?= $(TBB_DIR)/build/linux_armv7_gcc_cc4.9.2_libc2.19_kernel4.4.11_release
export LD_LIBRARY_PATH = $(TBB_DIR)/build/linux_armv7_gcc_cc4.9.2_libc2.19_kernel4.4.11_release

CPPFLAGS += -I $(TBB_INC_DIR)
LDFLAGS += -L $(TBB_LIB_DIR)


include ../makefile
