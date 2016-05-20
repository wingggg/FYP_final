SUBDIRS = LL optionGM utls extrema

all: 
	@for f in $(SUBDIRS) ; do $(MAKE) -C $$f all ; done
	ar crs libmser.a LL/LLconsis.o  LL/LLfile.o  LL/LLio.o  LL/LLmergeSort.o  LL/LL.o  LL/LLstr.o optionGM/optionDouble.o  optionGM/optionInt.o  optionGM/optionList.o  optionGM/optionMisc.o  optionGM/option.o utls/matrix.o  utls/timeutls.o extrema/boundary.o  extrema/getExtrema.o  extrema/optThresh.o   extrema/sortPixels.o  extrema/libExtrema.o  extrema/preprocess.o  extrema/suballoc.o

.PHONY: doc

doc:
	$(MAKE) -C doc doc

clean:
	@for f in $(SUBDIRS) ; do $(MAKE) -C $$f clean ; done

clean-doc:
	$(MAKE) -C doc clean

