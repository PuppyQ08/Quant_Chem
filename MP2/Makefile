EXENAME = MP2
OBJS = main.o MP2.o

CXX = clang++
CXXFLAGS = -std=c++1y -stdlib=libc++ -c -g -O2 -Wall -Wextra -pedantic
LD = clang++
LDFLAGS = -std=c++1y -stdlib=libc++ -lpthread -lm  -larmadillo

all : $(EXENAME)

$(EXENAME) : $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(EXENAME)

main.o : main.cpp scf.h
	$(CXX) $(CXXFLAGS) main.cpp

MP2.o : MP2.cpp MP2.h
	$(CXX) $(CXXFLAGS) MP2.cpp

clean :
	-rm -f *.o $(EXENAME)
