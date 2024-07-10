build:
	gcc -Wall -g -O2 -shared -fPIC -o libh5z_j2k.so -Lzstd/lib -Iwavelib/header -Lopenjpeg/build/bin -Iopenjpeg/src/lib/openjp2 -Iopenjpeg/build/src/lib/openjp2 -lzstd -lopenjp2 h5z_j2k.c wavelib/build/Bin/libwavelib.a

ubuntu:
	h5cc -Wall -g -O2 -shared -fPIC -o libh5z_j2k.so -Lzstd/lib -Iwavelib/header -Lopenjpeg/build/bin -Iopenjpeg/src/lib/openjp2 -Iopenjpeg/build/src/lib/openjp2 -I/usr/include/hdf5/serial/ -lzstd -lopenjp2 -lhdf5_serial h5z_j2k.c wavelib/build/Bin/libwavelib.a

setup:
	git submodule init
	git submodule update
	cd zstd && make
	mkdir openjpeg/build
	cd openjpeg/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make
	mkdir wavelib/build
	cd wavelib/build && cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON .. && make
