LIBRARIES= tbb

all:
	for LIBRARY in $(LIBRARIES) ; do \
	  cd $$LIBRARY ; make ; cd .. ; \
	done


clean:
	for LIBRARY in $(LIBRARIES) ; do \
	  cd $$LIBRARY ; make clean ; cd .. ; \
	done

