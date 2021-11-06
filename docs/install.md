# Install

The following is the installation sequence used on a small networked machine running debian jessie on 2017/02/06

This is after the network filesystem has been set up at `$NETDIR`:  

	sudo apt-get install build-essential gfortran cmake cmake-curses-gui git libfftw3-dev

Downloading tarballs into `$NETDIR`

	mkdir $NETDIR/mpich
	mkdir $NETDIR/metis
	mkdir $NETDIR/yaml

Get mpich from source (wget mpich-3.2)

	tar -xzf mpich-3.2.tar.gz
	./configure --prefix=$NETDIR/mpich 2>&1 | tee c.txt
	make
	mkdir $NETDIR/mpich
	sudo make install

Add to path in .bashrc

	# MPICH
	export PATH=$NETDIR/mpich/bin:$PATH
	export CPLUS_INCLUDE_PATH=$NETDIR/mpich/include:$CPLUS_INCLUDE_PATH
	export LIBRARY_PATH=$NETDIR/mpich/lib:$LIBRARY_PATH
	export LD_LIBRARY_PATH=$NETDIR/mpich/lib:$LD_LIBRARY_PATH

Get metis from source (wget metis-5.1.0)

	tar -xzf metis-5.1.0.tar.gz
	
Modify `include/metis.h`:  
	`#define IDXTYPEWIDTH 64`  
	`#define REALTYPEWIDTH 64`

Install `metis`:

	make config prefix=$NETDIR/metis
	make
	sudo make install

Get parmetis from source(wget parmetis-4.0.3)

	tar -xzf parmetis-4.0.3.tar.gz

Modify `metis/include/metis.h`:  
	`#define IDXTYPEWIDTH 64`  
	`#define REALTYPEWIDTH 64`
	
Install `parmetis`:

	make config prefix=$NETDIR/parmetis
	make
	sudo make install

Add to path in .bashrc

	# (PAR)METIS
	export PATH=$NETDIR/metis/bin:$PATH
	export CPLUS_INCLUDE_PATH=$NETDIR/metis/include:$CPLUS_INCLUDE_PATH
	export LIBRARY_PATH=$NETDIR/metis/lib:$LIBRARY_PATH
	export LD_LIBRARY_PATH=$NETDIR/metis/lib:$LD_LIBRARY_PATH

Get Charm++ from source (wget charm-6.7.1)

	tar -xzf charm-6.7.1.tar.gz
	mv $NETDIR/src/charm-6.7.1 $NETDIR/charm
	cd $NETDIR/charm
	./build charm++ mpi-linux-x86_64 -j8 -O3

Optionally:

	./build charm++ mpi-linux-x86_64 --with-production -j8 -O3

Add to path in .bashrc

	# CHARM
	export CHARMDIR=$NETDIR/charm
	export PATH=$CHARMDIR/bin:$PATH

Get yaml-cpp from github (so there's no boost dependencies)

	git clone https://github.com/jbeder/yaml-cpp.git
	mkdir build
	cd build
	ccmake ../

CMake settings:  
`CMAKE_INSTALL_PREFIX $NETDIR/yaml`

	make
	make install

Add to path in .bashrc

	# YAML
	export CPLUS_INCLUDE_PATH=$NETDIR/yaml/include:$CPLUS_INCLUDE_PATH
	export LIBRARY_PATH=$NETDIR/yaml/lib:$LIBRARY_PATH
	export LD_LIBRARY_PATH=$NETDIR/yaml/lib:$LD_LIBRARY_PATH

Get genet from github (fywang)

	git clone https://github.com/fywang/genet.git
	make -j4

Run a program:

	./charmrun +p4 -f $NETDIR/hostfile ./genet polynet.yml build
	./charmrun +p4 -f $NETDIR/hostfile ./genet polynet.yml part
	./charmrun +p4 -f $NETDIR/hostfile ./genet polynet.yml order

Get stacs from github (fywang)

	git clone https://github.com/fywang/stacs.git
	make -j8

Run a program:

	cd polynet
	cp ../../genet/polynet/* .
	mmv "polynet.o.*" "polynet.#1"
	cd ../
	./charmrun +p4 -f $NETDIR/hostfile ./stacs polynet.yml

Getting stacs with yarp needs yarp from github (robotology)
http://www.yarp.it/install_yarp_linux.html

	sudo apt-get install git cmake cmake-curses-gui libeigen3-dev libace-dev libedit-dev
	sudo apt-get install qtbase5-dev qtdeclarative5-dev qtmultimedia5-dev \
		qml-module-qtquick2 qml-module-qtquick-window2 \
		qml-module-qtmultimedia qml-module-qtquick-dialogs \
		qml-module-qtquick-controls libqt5svg5
	git clone https://github.com/robotology/yarp.git
	cd yarp
	mkdir build
	cd build
	ccmake ../

CMake settings:  
`CMAKE_INSTALL_PREFIX $NETDIR/yarp`  
`CREATE_GUIS, set to ON`  
`CREATE_LIB_MATH, set to ON`  
`CMAKE_SKIP_INSTALL_RPATH, set to ON`

	make
	sudo make install
	sudo ldconfig

Add to path in .bashrc

	# YARP
	export YARP_DIR=$NETDIR/yarp
	export YARP_DATA_DIRS=$YARP_DIR/share/yarp
	export PATH=$YARP_DIR/bin:$PATH
	export CPLUS_INCLUDE_PATH=$YARP_DIR/include:$CPLUS_INCLUDE_PATH
	export LIBRARY_PATH=$YARP_DIR/lib:$LIBRARY_PATH
	export LD_LIBRARY_PATH=$YARP_DIR/lib:$LD_LIBRARY_PATH

Rebuild stacs with yarp

	make clean;make -j8 WITHYARP=1

Run a program:

Terminal 1: `yarp server`  
Terminal 2: `./charmrun +p4 -f $NETDIR/hostfile ./stacs polynet.yml`  
Terminal 3: `yarp rpc /stacs/rpc`
