CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -g -fsanitize=address

SRC_DIR = src
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OUTPUT = raytracer

all: $(OUTPUT)

$(OUTPUT): $(SRC_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(OUTPUT)
