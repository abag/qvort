default: run

run:
	make -C src 
double:
	make -C src double
	#@make linkx
clean:
	make -C src clean
cleann:
	make -C src cleann
pristine: cleann
	make -C doc pristine
	rm -rf data/*
linkx:
	@for file in src/*.x; \
	do [ -e "`basename $$file`" ] || ln -s $$file .; \
	done
