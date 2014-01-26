default: run
run:
	make -C src 
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
	rm -f *.mat matlab.out
	rm -f *.vtk
	rm -f *.vdf *.pov *.log
	rm -rf unorm_data usup_data
	rm -f STOP
linkx:
	@for file in src/*.x; \
	do [ -e "`basename $$file`" ] || ln -s $$file .; \
	done
