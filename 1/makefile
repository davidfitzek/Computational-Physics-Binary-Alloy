# 
# Constants
#
CC = gcc # c-compiler
CFLAGS = -O3 -W # -O3 gives optimization, -W gives warnings
LIBS = -lm #-lgsl #-lgslcblas # libraries

OBJECTS = main.o # object files to produce
PROGRAM = main # name of executable to be built

#
# Rule that every .o files depends upon on the .c version of the file. Build .o files
#

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)
# $@ points to all .o files (left of the :) and $< to the first object to the right of : (in this case the corresponding .c file)

#
# Build everything that is needed to make a complete build.
#

all: $(PROGRAM) 

#
# Compilation
#

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
# $@ points to ${PROGRAM} (left of the :) and $^ to $(OBJECTS) (right of the :)

#
# make clean
#
clean:
	rm -f *.o
	touch *.c

