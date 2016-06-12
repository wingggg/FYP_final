  
CC = g++

ifeq ($(HOSTTYPE),intel-mac)
        RHOME = /Volumes/work
        MACFLAGS =  -arch x86_64
        MACINCLUDES = -I/opt/local/include -I/usr/X11/include
        MACDEF = -DMAC_OS
else
        RHOME = /user/cvsspstf/ees1km/Dropbox/code
        MACINCLUDES = -I/usr/include 
endif


LIBDIR = -L$(RHOME)/lib
INCDIR = -I$(RHOME)/lib
OBJDIR = ./Debug
TARGET = $(FILES:=.ln)
INCLUDES = $(INCDIR) $(INCDIR)/ann/include $(INCDIR)/mser/LL $(INCDIR)/mser/utls $(INCDIR)/mser/imageutls $(INCDIR)/mser/optionGM $(INCDIR)/face $(MACINCLUDES) $(INCDIR)/opencv/include/opencv $(INCDIR)/opencv/include/opencv2 $(INCDIR)/opencv/include

LFLAGS = -DWITH_LIBPNG -DWITH_LIBJPEG #-DWITH_LIBTIFF 
CFLAGS =  -c -O3 -Wall -fPIC -Wno-deprecated $(MACFLAGS) $(LFLAGS) $(MACDEF)


.SUFFIXES = .c .cc .cpp 

#C_FILES = $(FILES:=.c)
CPP_FILES = $(FILES:=.cpp)
O_FILES = $(FILES:%=$(OBJDIR)/%.o)
SHAREDLIB = $(LIB:=.so)

$(OBJDIR)/%.o: %.cc   
	$(CC) $(CFLAGS) $< -o $@

$(OBJDIR)/%.o: %.cpp   
	$(CC) $(CFLAGS) $< -o $@ $(INCLUDES)

$(LIB):  $(O_FILES)
	 ar rv $@ $?
	 ranlib $@

targets: $(LIB) 


