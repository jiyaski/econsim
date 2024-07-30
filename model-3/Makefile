
## makefile variables are just text-replacement 

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include

# compiler 
CXX = g++

# compiler flags 
CXXFLAGS = -Wall -g -I/usr/include/eigen3 -I$(INC_DIR) -Wno-ignored-attributes

# linker flags 
# LDFLAGS = 

# libraries to link 
LIBS = -lboost_program_options -lboost_math_tr1

# executable file to build 
TARGET = econsim

# list of source/object files 
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))


# Rule to link the object files into the executable
$(TARGET): $(OBJS) 
	@$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)


# Rule to compile the precompiled header
$(INC_DIR)/pch.h.gch: $(INC_DIR)/pch.h
	@$(CXX) $(CXXFLAGS) -x c++-header $(INC_DIR)/pch.h -o $(INC_DIR)/pch.h.gch

# Rules to compile the source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/pch.h.gch
	@$(CXX) $(CXXFLAGS) -c $< -o $@


# Commands/macros 
.PHONY: clean run

clean:
	@rm -f $(TARGET) $(OBJS) $(INC_DIR)/pch.h.gch

run: $(TARGET)
	@./$(TARGET)

