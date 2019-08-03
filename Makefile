CC=g++
EXE=screwDis

# 
# Compiles source files: src/*.cpp, cannot use subdirectories in src
# Output files to be made: Obj/*.o
SRC=$(wildcard src/*.cpp)
OBJ=$(patsubst src/%.cpp, Obj/%.o, $(SRC))

.SUFFIXES: .cpp .o

all: CFLAGS+= -O3
all: $(SRC) $(EXE)

debug: CFLAGS+= -g -ggdb
debug: LDFLAGS+= -g
debug: $(SRC) $(EXE)

Obj/%.o: src/%.cpp
	@mkdir -p Obj
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJ)
	rm -rf $(EXE)            
	$(CC) -lm $(OBJ) -o $@

clean:
	@rm -rf Obj
	rm -rf $(EXE)

