.PHONY: all examples clean

all: examples

examples: max_double_array.o

max_double_array.o: examples/max_double_array.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./src \
		examples/max_double_array.cpp -o max_double_array.o

clean:
	rm -rf `find . -name *.o`
