CXX = g++
CXXFLAGS = -O3 -Wall -pedantic -Wextra -Weffc++ -std=c++17 -Wshadow 

md: md.cpp Parser.o Atom.o Space.o Outer.o Utils.o
	$(CXX) $(CXXFLAGS) $< Atom.o Space.o Parser.o Outer.o Utils.o -o$@

Atom.o: Atom.cpp Atom.hpp
	$(CXX) $(CXXFLAGS) $< -c -o$@

Space.o: Space.cpp Space.hpp
	$(CXX) $(CXXFLAGS) $< -c -o$@

Parser.o: Parser/Parser.cpp Parser/Parser.hpp
	$(CXX) $(CXXFLAGS) $< -c -o$@

Outer.o: Outer/Out.cpp Outer/Out.hpp
	$(CXX) $(CXXFLAGS) $< -c -o$@

Utils.o: Utils/Utils.cpp Utils/Utils.hpp
	$(CXX) $(CXXFLAGS) $< -c -o$@

removeVTK:
	rm -f vtk/*.vtk

clean:
	rm -f *.o md log*