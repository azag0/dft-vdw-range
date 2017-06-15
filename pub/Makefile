BLDDIR = _build

default: revtex

revtex: $(BLDDIR)/preprint-revtex.pdf

preprint: $(BLDDIR)/preprint.pdf

$(BLDDIR)/%.pdf: %.tex FORCE
	env max_print_line=1000000 latexmk -pdf -outdir=$(BLDDIR) -interaction=nonstopmode $<

bib: refs.bib

refs.bib: FORCE
	-curl -f -o $@ "http://localhost:23119/better-bibtex/collection?/0/References/mbd-scan.bibtex&exportNotes=false&useJournalAbbreviation=true"

FORCE:
