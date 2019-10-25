CC   := g++
INC  := -I include
LIB  := 

Program := Liscio
OBJ  := main.o base.o func.o easy.o nve.o

VPATH = src lib/prime lib/math lib/hes lib/md

.PHONY: all build_msg clean

all: $(Program)

build_msg:
	@printf "#\n# Building $(Program)\n#\n"

clean:
	rm -f *.o

$(OBJ): %.o : %.cpp
	$(CC) -c $< -o $@ $(INC)

$(Program): build_msg $(OBJ)
	$(CC) $(OBJ) -o $@ $(INC) $(LIB)
	ls -hl $(Program)
	size $(Program)
