
CC = g++ -std=c++11 -pthread

ifeq ($(HOSTTYPE),intel-mac)
        RHOME = /Users/kma/Dropbox/code
        MACFLAGS =  -arch x86_64 
        MACINCLUDES = -I/opt/local/include 
        #MACINCLUDES = -I/usr/X11/include 
        MACDEF = -DMAC_OS
		BIN_EXT = ln 
else
        RHOME = /home/pi/Documents/FYP_final/data_pi/bcode
        #/home/cvssp/stf/ees1km/Dropbox/code
        #RHOME = /home/kma/Dropbox/code#/user/cvsspstf/ees1km/Dropbox/code
        BIN_EXT = ln
endif




LIBDIR = -L$(RHOME)/lib
INCDIR = -I$(RHOME)/lib
OBJDIR = ./Debug
TARGET = $(FILES:=.$(BIN_EXT))
#INCLUDES = $(INCDIR)/ann/include $(INCDIR)/mser/LL $(INCDIR)/mser/utls $(INCDIR)/mser/imageutls $(INCDIR)/mser/optionGM $(INCDIR)/face $(MACINCLUDES)  $(INCDIR)/opencv/include/opencv $(INCDIR)/opencv/include -I/usr/local/include -I/usr/include

LFLAGS = -DWITH_LIBPNG -DWITH_LIBJPEG #-DWITH_LIBTIFF 
#CFLAGS =  -c -O3 -fopenmp -Wall -fPIC -Wno-deprecated $(MACFLAGS) $(LFLAGS) $(MACDEF)
CFLAGS =  -c -O3 -Wall -fPIC -Wno-deprecated $(LDFLAGS) $(LFLAGS) $(MACDEF)

#LIBRARY += -L/opt/local/lib -lpng  -ltiff  -ljpeg -lz -lm  -fopenmp
LIBRARY +=  -lpng -ljpeg -lz -lm 


vpath %.o $(OBJDIR)

HOME_BIN=  $(RHOME)/bin

C_FILES = $(FILES:=.c)
O_FILES = $(FILES:%=$(OBJDIR)/%.o)


$(OBJDIR)/%.o: %.cc  
	$(CC) $(CFLAGS) $< -o $@

$(OBJDIR)/%.o: %.cpp $(DEPEND)
	$(CC) $(CFLAGS) $< -o $@  $(INCLUDES) -pthread -std=c++11

$(TARGET):  $(O_FILES)
	 $(CC)  -o $(HOME_BIN)/$@ $(O_FILES) $(LDFLAGS) -lpthread -pthread -std=c++11 $(LIBRARY) -lpthread -pthread -std=c++11
#	 $(CC)  -o $(HOME_BIN)/$@ $(O_FILES) $(MACFLAGS) -lpthread $(LIBRARY) -fopenmp -lpthread

# scp $(HOME_BIN)/$@ ees1km@131.227.95.122:/vol/vssp/recog/bin/.


all:  $(TARGET)
targets: $(TARGET)

clean: 
	rm -fr Debug/*
