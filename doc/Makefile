
# get ms2 version from file, store in environment variable for use in Dixyfile
export MS2_VERSION := $(shell grep -oP "'\K[^']+" ../src/ms2_version.F90)

all:
	@make clean
	@make html
	cd latex; make pdf;  cd ..

html: clean
	@doxygen ./Doxyfile

latex: clean
	@doxygen ./Doxyfile
	cd latex; make pdf;  cd ..

clean:
	rm -rf ./html ./latex
