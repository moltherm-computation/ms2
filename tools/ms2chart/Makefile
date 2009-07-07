.PHONY: all clean clean_javadoc_dir clean_class_dir

all: ms2chart.jar ms2chart-javadocs.tar.bz2

ms2chart.jar: ms2package manifest.mf
	jar cmf manifest.mf ms2chart.jar ms2package/*.class org
	chmod a+r ms2chart.jar

ms2package: Ms2chart.java ShowFileFrame.java ShowTableFrame.java
	rm -fr ms2package
#	javac -d . *.java
	javac -classpath .:lib/swing-layout-1.0.3.jar:lib/libs.zip -d . *.java

ms2chart-javadocs.tar.bz2: ms2chart-javadocs
	tar -cvjf ms2chart-javadocs.tar.bz2 ms2chart-javadocs
	chmod a+r ms2chart-javadocs.tar.bz2
	
ms2chart-javadocs: Ms2chart.java ShowFileFrame.java ShowTableFrame.java
	rm -fr ms2chart-javadocs
	mkdir -p ms2chart-javadocs
	javadoc -d ms2chart-javadocs -version -author Ms2chart.java ShowFileFrame.java ShowTableFrame.java
	chmod -R a+rX ms2chart-javadocs

clean: clean_javadoc_dir clean_class_dir
	rm -fr ms2chart-javadocs.tar.bz2 ms2chart.jar ms2package ms2chart-javadocs

clean_javadoc_dir:
	rm -fr ms2chart-javadocs

clean_class_dir:
	rm -fr ms2package
