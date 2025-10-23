CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Iinclude
LDFLAGS =

SRC_DIR = src
INC_DIR = include
EXAMPLES_DIR = examples
TEST_DIR = tests
BIN_DIR = bin

SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)

EXAMPLES = pendulum_simulation harmonic_oscillator planetary_motion heat_equation
EXAMPLE_BINS = $(addprefix $(BIN_DIR)/, $(EXAMPLES))

.PHONY: all clean examples tests

all: $(BIN_DIR) $(OBJECTS) examples tests

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

examples: $(EXAMPLE_BINS)

$(BIN_DIR)/%: $(EXAMPLES_DIR)/%.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) -o $@ $(LDFLAGS)

tests: $(BIN_DIR)/unit_tests
	./$(BIN_DIR)/unit_tests

$(BIN_DIR)/unit_tests: $(TEST_DIR)/unit_tests.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -rf $(BIN_DIR)/*.o $(BIN_DIR)/$(EXAMPLES) $(BIN_DIR)/unit_tests
	rm -f *.csv *.png

run-all: examples
	@echo "Ejecutando todas las simulaciones..."
	@./$(BIN_DIR)/pendulum_simulation
	@./$(BIN_DIR)/harmonic_oscillator
	@./$(BIN_DIR)/planetary_motion
	@./$(BIN_DIR)/heat_equation

