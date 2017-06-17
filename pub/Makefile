BLDDIR = _build
export max_print_line = 1000000

default: revtex

revtex: $(BLDDIR)/preprint-revtex.pdf

preprint: $(BLDDIR)/preprint.pdf

all: revtex preprint figs

figs:
	@mkdir -p $(BLDDIR)
	pdflatex -output-directory=$(BLDDIR) -shell-escape export-figures.tex

$(BLDDIR)/%.pdf: %.tex refs-zotero.bib FORCE
	latexmk -pdf -outdir=$(BLDDIR) -interaction=nonstopmode $<

bib: refs-zotero.bib

refs-zotero.bib: FORCE
	-curl -sfo $@ "http://localhost:23119/better-bibtex/collection?/0/References/mbd-scan.bibtex&exportNotes=false&useJournalAbbreviation=true"

clean:
	rm -rf $(BLDDIR)

FORCE:
