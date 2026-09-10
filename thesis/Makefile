# make        -> build main.pdf
# make clean  -> remove build artifacts
main.pdf: main.tex refs.bib
	latexmk -pdf -interaction=nonstopmode main.tex

clean:
	latexmk -C
	rm -f *.bbl *.blg
.PHONY: clean
