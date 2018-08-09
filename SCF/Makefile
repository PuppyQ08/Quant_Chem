EXENAME = scf
OBJS = main.o scf.o

CXX = clang++
CXXFLAGS = -std=c++1y -stdlib=libc++ -c -g -O2 -Wall -Wextra -pedantic
LD = clang++
LDFLAGS = -std=c++1y -stdlib=libc++ -lpthread -lm  -larmadillo

all : $(EXENAME)

$(EXENAME) : $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(EXENAME)

main.o : main.cpp scf.h
	$(CXX) $(CXXFLAGS) main.cpp

scf.o : scf.cpp scf.h
	$(CXX) $(CXXFLAGS) scf.cpp

clean :
	-rm -f *.o $(EXENAME)
