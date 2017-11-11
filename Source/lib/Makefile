CXX=g++4.8
LIBRARIES= tbb seq-align bwa SHD edlib opal/src Complete-Striped-Smith-Waterman-Library/src bwa-0.5.9

#edlib requires cmake.

all:
	for LIBRARY in $(LIBRARIES) ; do \
	  cd $$LIBRARY ; make ; cd .. ; \
	done


clean:
	for LIBRARY in $(LIBRARIES) ; do \
	  cd $$LIBRARY ; make clean ; cd .. ; \
	done

