include ${PETSC_DIR}/lib/petsc/conf/variables

CC=g++
VPATH=src
CFLAGS=-c -g -O2 -std=c++11 -I./include -I${PETSC_DIR}/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include
OBJ := $(patsubst $(VPATH)/%.cpp,%.o,$(wildcard $(VPATH)/*.cpp))
CPP := $(patsubst $(VPATH)/%.cpp,%.cpp,$(wildcard $(VPATH)/*.cpp))

TARGET=main

$(TARGET): $(OBJ)
	$(CC) $^ -o $@ $(PETSC_LIB)
	rm *.o

$(OBJ): $(CPP)
	$(CC) $(CFLAGS) $^

clean:
	rm -f $(TARGET)
