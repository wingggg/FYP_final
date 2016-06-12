include ../config
	
INCLUDE_PATHS = -I../LL -I../utls -I../imageutls -I../optionGM -I/usr/include/objc
CPP_FILES  = libExtrema.cpp sortPixels.cpp getExtrema.cpp optThresh.cpp boundary.cpp \
	  			 preprocess.cpp suballoc.cpp

OBJ_FILES = $(CPP_FILES:.cpp=.o)
LIB_NAME = libExtrema.a

all: $(LIB_NAME)

extrema.o: CXXFLAGS += $(IMAGE_LIBS_CFLAGS)

extrema-bin: $(LIB_NAME) extrema.o
	$(CXX) $(CXXFLAGS) extrema.o -L. -L../LL -L../utls -L../optionGM -lExtrema -lLL -lutls -lOptionGM $(IMAGE_LIBS_LFLAGS) -lrt -lz -lm -lpthread -static -o extrema-bin

$(LIB_NAME): $(OBJ_FILES)
	 $(AR) crs $(LIB_NAME) $(OBJ_FILES)

dep:
	makedepend -Y -f- -- $(CFLAGS) $(CXXFLAGS) $(INCLUDE_PATHS) -- $(CPP_FILES) extrema.cpp 2> /dev/null > depend

-include depend

GARBAGE = $(OBJ_FILES) extrema-bin extrema.o core *~ *.mexlx *.mexglx *.a
clean:
	$(RM) $(GARBAGE)
