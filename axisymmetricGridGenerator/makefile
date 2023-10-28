# define variables
HDRDIR  = ./include

# specify compilers to use

CC = g++
LD = g++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I$(HDRDIR)  -O3 -std=c++11 -g3

# link flags to be used
LDFLAGS	= -I$(HDRDIR) 

# libraries to be linked in
LIBS = -lm

# types of files we are going to construct rules for
.SUFFIXES: .cpp

# rule for .cpp files
.cpp.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.cpp

# list of objects to be compiled
OBJS = \
src/main.o\
src/includedSourceFiles.o\
#src/XXXXXXXXXXXXXXXXXXXXXX.o\
#src/XXXXXXXXXXXXXXXXXXXXXX.o\


all:$(OBJS)
	$(LD)  $(LDFLAGS) -o main $(OBJS) $(LIBS)
	rm -r $(OBJS)

# what to do if user types "make clean"
clean :
	rm main

realclean :
	rm -r $(OBJS)
