PKGS=sdl2 glew gl libpng16
CC=clang
#OPT=-Ofast
OPT=-O0 -ggdb3
CCOMMON=$(OPT) -Wall $(shell pkg-config $(PKGS) --cflags)
CFLAGS=--std=c99 $(CCOMMON)
LINK=-lm $(shell pkg-config $(PKGS) --libs)
GLSL2INC=./glsl2inc.pl

all: main

a.o: a.c
	$(CC) $(CFLAGS) -c a.c

m.o: m.c
	$(CC) $(CFLAGS) -c m.c

shader0.glsl.inc: shader0.vert.glsl shader0.frag.glsl
	$(GLSL2INC) shader0 shader0.glsl.inc shader0.vert.glsl shader0.frag.glsl

main.o: main.c shader0.glsl.inc
	$(CC) $(CFLAGS) -c main.c

main: main.o a.o m.o
	$(CC) main.o a.o m.o -o main $(LINK)

clean:
	rm -rf *.o *.glsl.inc main

