CC   := g++
INC  := -I include -I ~/Documents/Development/usr/include
LIB  := 

EXE  := Liscio
OBJ  := main.o base.o func.o easy.o nve.o rho.o

VPATH = src lib/prime lib/math lib/hes lib/md lib/density

.PHONY: all build_msg install clean

all: $(EXE)

build_msg:
	@printf "#\n# Building $(EXE)\n#\n"

$(OBJ): %.o : %.cpp
	$(CC) -c $< -o $@ $(INC) $(LIB)

$(EXE): build_msg $(OBJ)
	$(CC) $(OBJ) -o $@ $(INC) $(LIB)
	ls -hl $(EXE)
	size $(EXE)

install:
	cp $(EXE) ~/bin/$(EXE)
	@printf "Normal installation. Congratulations!\n"

clean:
	rm -f *.o