PROJNAME := note

# Uncomment these lines to compile outside of tara
#export TEXINPUTS:=$(PWD)/../template:.:
#export BSTINPUTS:=$(PWD)/../template:.:
#export BIBINPUTS:=$(PWD)/../bib


main:
	pdflatex $(PROJNAME)

bib:
	bibtex $(PROJNAME)

clean:
	rm -f *.aux *.out *.bbl *.blg *.toc *.lof *.log *.ptp *.pdf *.sty

