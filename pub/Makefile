BLDDIR = build

all: preprint suppl

preprint: refs.bib | blddir
	latexmk -pdf -halt-on-error -jobname=$(BLDDIR)/$@ $@
	rsync -a $(BLDDIR)/$@.pdf .

suppl: | blddir
	latexmk -pdf -halt-on-error -jobname=$(BLDDIR)/$@ $@
	rsync -a $(BLDDIR)/$@.pdf .

bib: refs.bib

refs.bib: FORCE
	-curl -f -o $@ "http://localhost:23119/better-bibtex/collection?/0/References/mbd-scan.bibtex&exportNotes=false&useJournalAbbreviation=true"

FORCE:

blddir:
	mkdir -p $(BLDDIR)
	@ln -fns . $(BLDDIR)/$(BLDDIR)

