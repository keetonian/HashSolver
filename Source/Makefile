FOLDERS=L1Baseline L2Baseline

all: 
	for FOLDER in $(FOLDERS) ; do \
	  cd $$FOLDER ; make ; cd .. ; \
	done ;\
	cd lib; make; cd .. ;

clean:
	for FOLDER in $(FOLDERS) ; do \
		cd $$FOLDER ; make clean ; cd .. ; \
	done
