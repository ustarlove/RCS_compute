CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall
LDFLAGS = -lm

TARGET = sphere_rcs
SRCS = main.cpp mesh.cpp rwg.cpp efie.cpp cg_solver.cpp farfield.cpp \
       fmm.cpp fmm_translator.cpp mlfma.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run