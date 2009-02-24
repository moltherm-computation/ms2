.PHONY: all clean clean_class_dir

all: ms2par.jar

ms2par.jar: inputfilegen manifest.txt
	jar cmf manifest.txt ms2par.jar inputfilegen/*.class lib/swing-layout-1.0.jar
	chmod a+r ms2par.jar

inputfilegen: src/inputfilegen/*.java
	rm -fr inputfilegen
	javac -classpath lib/swing-layout-1.0.jar -d . src/inputfilegen/*.java



clean: clean_class_dir
	rm -fr ms2par.jar inputfilegen

clean_class_dir:
	rm -fr inputfilegen
