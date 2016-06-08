include ../config

INCLUDE_PATHS = -I.

LIB = libLL.a 
C_FILES = LL.c LLio.c LLmergeSort.c LLconsis.c LLstr.c LLfile.c
OBJ_FILES = $(C_FILES:.c=.o)

all: $(LIB)

$(LIB): $(OBJ_FILES)
	$(AR) crs $(LIB) $(OBJ_FILES)

GARBAGE = $(OBJ_FILES) $(LIB) core *~

.PHONY:clean
clean:
	$(RM) $(GARBAGE)
