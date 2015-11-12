include makefile.in

HDR = $(wildcard *.h)
SRC = $(wildcard *.cc)
SRC += ../extern/fparser/fparser.cc
OBJ = $(patsubst %.cc,%.o,$(SRC))
INC = -I../extern/fparser


TARGET = aim_run_extended

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ) $(LDFLAGS)

$(OBJ): $(HDR)

.cc.o: $(HDR)
	$(CXX) $(CFLAGS) $(INC) -c $*.cc -o $*.o

clean:
	rm -f $(OBJ) $(TARGET)

allclean:
	make clean
