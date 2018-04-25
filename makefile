SRCDIR= ./src
install:
	cd $(SRCDIR); make

macos_install:
	cd $(SRCDIR); make -f makefile.macos

boost:
	cd boost_1_65_0; chmod +x install_boost.sh; ./install_boost.sh
	
clean:
	cd $(SRCDIR); make clean


