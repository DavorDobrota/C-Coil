install:
	CPPFLAGS="-std=c++17 -O3" python3 setup.py install

build:
	CPPFLAGS="-std=c++17 -O3" python3 setup.py build

clean:
	rm -rf build/
