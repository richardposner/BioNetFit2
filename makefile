SRCDIR= ./src
install:
	cd $(SRCDIR); make

boost:
	cd boost_1.59.0; chmod +x install_boost.sh; ./install_boost.sh
	
clean:
	cd $(SRCDIR); make clean


