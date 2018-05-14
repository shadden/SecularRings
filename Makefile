CC=gcc
INCLUDE=./
CFLAGS= -std=c99
LFLAGS=  -lm -lgsl -lgslcblas 
%.o: %.c 
	gcc -c $< $(CFLAGS)
testme: testme.c orbits.o vectors.o averaging.o
	$(CC) -I$(INCLUDE) $^ -o $@ $(CFLAGS) $(LFLAGS)
	./$@
testme2: testme2.c orbits.o vectors.o averaging.o
	$(CC) -I$(INCLUDE) $^ -o $@ $(CFLAGS) $(LFLAGS)
	./$@
clean:
	rm -f *.o *.mod 

