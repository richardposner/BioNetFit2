BINDIR= ./bin
#Razi added to compile everytime
#Message:
#	@echo 'running main makefile'

#cleanup:
#	rm -f ./bin/*.exe
#	rm -f ./bin/src/*.d
#	rm -f ./bin/src/*.o
#	rm -f ./bin/src/code/*.o
#	rm -f ./bin/src/code/*.d

install:
	@echo 'running /bin/main makefile'
	cd $(BINDIR); make

#razi commented, only one run per pc is needed
boost:
	cd boost_1.59.0; chmod +x install_boost.sh; ./install_boost.sh
	
clean:
	cd $(BINDIR); make clean


