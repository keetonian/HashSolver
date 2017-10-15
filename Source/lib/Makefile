LIBRARIES= tbb seq-align bwa Shifted-Hamming-Distance edlib

all:
	for LIBRARY in $(LIBRARIES) ; do \
	  cd $$LIBRARY ; make ; cd .. ; \
	done


clean:
	for LIBRARY in $(LIBRARIES) ; do \
	  cd $$LIBRARY ; make clean ; cd .. ; \
	done

