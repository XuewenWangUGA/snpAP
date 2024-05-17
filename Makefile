

ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

all: 
	@cd Bam2snpAD && make
	@cd snpAD && make

clean: 
	@cd Bam2snpAD && make clean
	@cd snpAD && make clean

install: all
	install -d $(DESTDIR)$(PREFIX)/bin/
	install -m 755 Bam2snpAD/Bam2snpAD $(DESTDIR)$(PREFIX)/bin/
	install -m 755 Bam2snpAD/snpADjoin $(DESTDIR)$(PREFIX)/bin/
	install -m 755 snpAD/snpAD $(DESTDIR)$(PREFIX)/bin/
	install -m 755 snpAD/snpADCall $(DESTDIR)$(PREFIX)/bin/
	
release:
	echo "#define SNPAD_VERSION \"`cat VERSION`\"" > version.h
	git commit version.h -m "new version"
	git tag -a `cat VERSION` -m "release `cat VERSION`"
	git archive --format=tar.gz --prefix=snpAD-`git describe`/ `git describe` >../snpAD-`git describe`.tar.gz
