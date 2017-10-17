LIBRARIES= tbb seq-align bwa Shifted-Hamming-Distance edlib opal/src Complete-Striped-Smith-Waterman-Library/src

#edlib requires cmake.

all:
	for LIBRARY in $(LIBRARIES) ; do \
	  cd $$LIBRARY ; make ; cd .. ; \
	done


clean:
	for LIBRARY in $(LIBRARIES) ; do \
	  cd $$LIBRARY ; make clean ; cd .. ; \
	done

