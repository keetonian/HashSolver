FOLDERS=lib L1Baseline L2Baseline

all: 
	for FOLDER in $(FOLDERS) ; do \
	  cd $$FOLDER ; make ; cd .. ; \
	done

clean:
	for FOLDER in $(FOLDERS) ; do \
		cd $$FOLDER ; make clean ; cd .. ; \
	done
