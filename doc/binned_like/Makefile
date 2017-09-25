manuscript=binned_like

TEXFILES = $(wildcard *.tex $(manuscript).bib )

all: $(manuscript).pdf

plots:
	python scripts/efficiency.py
	python scripts/countmap.py

pdf: $(manuscript).pdflatex

ps: $(manuscript).ps

ps2pdf: $(manuscript).ps
	pstopdf $(manuscript).ps

force: touch all

bib: 
	bibtex $(manuscript)

%.pdf: %.dvi
	dvipdf $<

%.ps: %.dvi
	dvips $<

%.dvi: %.tex $(TEXFILES)
	latex $(manuscript).tex
	#bibtex $(manuscript)
	#latex $(manuscript).tex
	latex $(manuscript).tex

%.pdflatex: %.tex $(TEXFILES)
	pdflatex $(manuscript).tex
	#bibtex $(manuscript)
	#pdflatex $(manuscript).tex
	pdflatex $(manuscript).tex

touch:
	touch $(manuscript).tex

clean:
	rm -f *.out *.log *.dvi *.aux *~ */*~ *.blg *.bbl *.bb *.nav $(manuscript).pdf *.snm *.toc $(manuscript).ps 

cleanall:
	make clean
	rm -f *.pdf
