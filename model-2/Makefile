.PHONY: clean run

# Variables are just text-replacement 
SRCDIR = src
CXX = g++
CXXFLAGS = -Wall -g
TARGET = econsim
OBJS = $(SRCDIR)/main.o $(SRCDIR)/distributions.o

# Rule to link the object files into the executable
$(TARGET): $(OBJS)
	@$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule to compile the source files into object files
$(SRCDIR)/main.o: $(SRCDIR)/main.cpp $(SRCDIR)/distributions.h
	@$(CXX) $(CXXFLAGS) -c $(SRCDIR)/main.cpp -o $(SRCDIR)/main.o

$(SRCDIR)/distributions.o: $(SRCDIR)/distributions.cpp $(SRCDIR)/distributions.h
	@$(CXX) $(CXXFLAGS) -c $(SRCDIR)/distributions.cpp -o $(SRCDIR)/distributions.o

# Commands/macros 
clean:
	@rm -f $(TARGET) $(OBJS)

run: $(TARGET)
	@./$(TARGET)

