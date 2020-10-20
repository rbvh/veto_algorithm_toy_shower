CXXFLAGS= -std=c++11 -O -I. `root-config --cflags --glibs`

All:ShowerRunner

ShowerRunner: ShowerRunner.cc ShowerBase.h ShowerVetoMaxDumb.h ShowerVetoMax.h ShowerMaxVetoDumb.h ShowerMaxVeto.h ShowerGenerateSelectDumb.h ShowerGenerateSelect.h ShowerMaxVetoSmart.h Basics.o PythiaStdlib.h
	$(CXX) $(CXXFLAGS) -o $@ ShowerRunner.cc Basics.o

.PHONY: clean

clean:
	rm -f ShowerRunner
	rm -f Basics.o

