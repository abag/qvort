

LATEX=latex
GNUPLOT=gnuplot

default=linedensity

linedensity:linedensity.gp
	make clean
	$(GNUPLOT) linedensity.gp
	$(LATEX) figure
	$(LATEX) figure
	$(LATEX) figure
	dvips figure.dvi
	ps2eps figure.ps
	rm figure.dvi figure.log figure.ps
	mv figure.eps linedensity.eps
	epstopdf linedensity.eps
curvature:curvature.gp
	make clean
	$(GNUPLOT) curvature.gp
	$(LATEX) figure
	$(LATEX) figure
	$(LATEX) figure
	dvips figure.dvi
	ps2eps figure.ps
	rm figure.dvi figure.log figure.ps
	mv figure.eps curvature.eps
	epstopdf curvature.eps
clean:
	rm -f figure.eps figure.ps figure.dvi figure.aux figure.tex
pristine:
	make clean
	rm -f *.eps *.pdf 


