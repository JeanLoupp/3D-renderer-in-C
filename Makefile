FILEPATH = "models/circuit_mario/circuit_tri.obj"
WIDTH = 1200
HEIGHT = 800
FOCAL = 1000
DISTANCE = 25

SIMPLE = 1

compile:
	gcc main.c -o main $(shell sdl2-config --cflags --libs) -lSDL2_image -O2
	./main $(FILEPATH) $(WIDTH) $(HEIGHT) $(FOCAL) $(DISTANCE)

run:
	./main $(FILEPATH) $(WIDTH) $(HEIGHT) $(FOCAL) $(DISTANCE)