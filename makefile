SRCDIR= ./src
install:
	cd $(SRCDIR); make

macos_install:
	cd $(SRCDIR); make -f makefile.macos

clang_install:
	cd $(SRCDIR); make -f makefile.clang

ubuntu_install:
	cd $(SRCDIR); make -f makefile.ubuntu

c17_install:
	cd $(SRCDIR); make -f makefile.c17

boost:
	cd boost_1_65_0; chmod +x install_boost.sh; ./install_boost.sh

clean:
	cd $(SRCDIR); make clean

docker_clean:
	cd $(SRCDIR); make docker_clean


