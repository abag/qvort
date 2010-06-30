#$Id: Makefile 58 2010-02-24 23:04:25Z abag $
LATEX=pdflatex
BIBTEX=bibtex
VIEWER=acroread

default: manual

manual: manual.tex Makefile
	$(LATEX) manual
# 	$(BIBTEX) manual
#	$(LATEX) manual
	$(LATEX) manual
	$(VIEWER) manual.pdf
clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.blg *~
pristine: clean
	rm -f *.pdf *.dvi *.bbl
    

