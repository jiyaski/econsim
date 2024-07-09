.PHONY: clean run

## makefile variables are just text-replacement 

SRCDIR = src

# compiler 
CXX = g++ 

# compiler flags 
CXXFLAGS = -Wall -g -I/usr/include/eigen3

# linker flags 
# LDFLAGS = 

# libraries to link 
LIBS = -lboost_program_options 

# executable file to build 
TARGET = econsim

# list of object files 
OBJS = $(SRCDIR)/main.o 


# Rule to link the object files into the executable
$(TARGET): $(OBJS)
	@$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)


# Rule to compile the source files into object files
$(SRCDIR)/main.o: $(SRCDIR)/main.cpp 
	@$(CXX) $(CXXFLAGS) -c $(SRCDIR)/main.cpp -o $(SRCDIR)/main.o


# Commands/macros 
clean:
	@rm -f $(TARGET) $(OBJS)

run: $(TARGET)
	@./$(TARGET)

