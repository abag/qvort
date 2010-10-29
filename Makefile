default: run

run:
	make -C src 
double:
	make -C src double
man:
	make -C doc
clean:
	make -C src clean
cleann:
	make -C src cleann
pristine: cleann
	make -C doc pristine
	rm -rf data/*
	rm -f *.eps *.png
	rm -f *.mat
linkx:
	@for file in src/*.x; \
	do [ -e "`basename $$file`" ] || ln -s $$file .; \
	done
