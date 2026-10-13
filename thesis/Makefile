# make           -> build main.pdf (full thesis) and intro_methods.pdf (partial submission)
# make clean     -> remove build artifacts
SECTIONS := $(wildcard sections/*.tex) preamble.tex refs.bib

all: main.pdf intro_methods.pdf

main.pdf: main.tex $(SECTIONS)
	latexmk -pdf -interaction=nonstopmode main.tex

intro_methods.pdf: intro_methods.tex $(SECTIONS)
	latexmk -pdf -interaction=nonstopmode intro_methods.tex

clean:
	latexmk -C main.tex intro_methods.tex
	rm -f *.bbl *.blg
.PHONY: all clean
