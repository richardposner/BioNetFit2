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
	cd $(SRCDIR); OPTIONS=$(OPTIONS) make -f makefile.c17

monsoon_install:
	cd $(SRCDIR); make -f makefile.monsoon

boost_1_67_0.tmp.downloaded:
	curl -L https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.gz > /tmp/boost_1_67_0.tar.gz && touch boost_1_67_0.tmp.downloaded

boost_1_67_0.tar.gz:
	curl -L https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.gz > boost_1_67_0.tar.gz

boost_1_67_0.expanded: boost_1_67_0.tar.gz
	tar xzf boost_1_67_0.tar.gz && touch boost_1_67_0.expanded

boost_1_67_0.configured: boost_1_67_0.expanded
	cd boost_1_67_0 && bash ./bootstrap.sh --prefix=../ --with-libraries=filesystem,iostreams,mpi,regex,system,program_options,serialization && echo "using mpi ;" >> project-config.jam && cd - && touch boost_1_67_0.configured

boost: boost_1_67_0.configured
	cd boost_1_67_0 && ./b2 install

docker_boost: boost_1_67_0.tmp.downloaded
	cd /tmp && tar xzf boost_1_67_0.tar.gz && cd - && cp -r /tmp/boost_1_67_0 . && cd boost_1_67_0 && bash ./bootstrap.sh --prefix=../ --with-libraries=filesystem,iostreams,mpi,regex,system,program_options,serialization && echo "using mpi ;" >> project-config.jam && ./b2 install

clean:
	cd $(SRCDIR); make clean

docker_clean:
	cd $(SRCDIR); make docker_clean

docker_boost_clean:
	pwd; cd .; pwd; rm -rvf boost_1_67_0; rm -rf /tmp/boost_1_67_0*

