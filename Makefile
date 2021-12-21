compiler = icpc
flags = 
flags = -I. -I./eigen/ -lm -std=c++17

headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)

executables: solver

%.o: %.cpp $(headers)
	$(compiler) -c -o $@ $< $(flags)

solver: $(objects)
	$(compiler) -o $@ $^ $(flags)

clean:
	rm -f *.o

run:
	./solver