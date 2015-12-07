
CXXFLAGS = -Wall -Wextra -Wno-unused-parameter -O3
#-Wno-write-strings
LDLIBS = -lglsc3d -lGL -lglut -lGLU -lpng -I/usr/include -I/usr/X11R6/include /usr/lib/libglsc3d.a -L/opt/X11/lib

ALL = main
RUN = ./main

all : ${ALL}

clean :
	rm -rf *.o ${ALL}

run : ${ALL}
	${RUN}

main : vertex.o cellToVert.o myglsc.o vertToVert.o update.o vertToCell.o readfile.o makePolygon.o writeMesh.o functions.o
