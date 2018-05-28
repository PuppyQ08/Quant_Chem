EXENAME = molec_geo
OBJS = main.o molec_geo.o

CXX = clang++
CXXFLAGS = -std=c++1y -stdlib=libc++ -c -g -O0 -Wall -Wextra -pedantic
LD = clang++
LDFLAGS = -std=c++1y -stdlib=libc++ -lpthread -lm

all : $(EXENAME)

$(EXENAME) : $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(EXENAME)

main.o : main.cpp molec_geo.h
	$(CXX) $(CXXFLAGS) main.cpp

molec_geo.o : molec_geo.cpp molec_geo.h
	$(CXX) $(CXXFLAGS) molec_geo.cpp

clean :
	-rm -f *.o $(EXENAME)
