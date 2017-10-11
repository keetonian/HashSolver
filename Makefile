all:
	cd Source/lib && make && cd .. && make && cd ../..; \
	cd bin; \
	sh ./link_binaries.sh; \
	cd ..;

clean:
	cd Source && make clean && cd ..; \
	find ./bin/ -type l -exec rm {} \;
