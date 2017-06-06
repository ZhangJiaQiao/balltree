CC := g++
FLAGS := -std=c++11 -w
SRC_DIR := src
BUILD_DIR := build
BIN_DIR := bin
INCLUDE_DIR := include
INCLUDE := -I./$(INCLUDE_DIR)

$(BIN_DIR)/main: $(BUILD_DIR)/test.o $(BUILD_DIR)/BallTree.o $(BUILD_DIR)/Utility.o
	@mkdir -p $(BIN_DIR)
	$(CC) $(FLAGS) $(INCLUDE) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) $(FLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@rm -rf $(BUILD_DIR)
	@rm -rf $(BIN_DIR)