
CC = gcc
CFLAGS = -O3 -Wall # -O3 # -O0 for valgrind
DBGFLAGS = 
LIBS = -lm -lgsl -lgslcblas

HEADERS = 
OBJECTS = main.o
PROGRAM = main


%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS) $(DBGFLAGS)

all: $(PROGRAM)

debug: DBGFLAGS += -DDEBUG -g -O0
debug: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(DBGFLAGS)

clean:
	rm -f *.o
	touch *.c

