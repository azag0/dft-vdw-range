FIGS = range-curves.pdf s66-dists.pdf mbd-param-fitting.pdf \
	   vv10-param-fitting.pdf d3-param-fitting.pdf bz2.png \
	   cor21.png cor22.png gr2.png flakes.pdf scan-interp.pdf
STANDALONE = $(patsubst %,standalone-figures/fig-%.pdf,1 2 3 4 5)
TOCOPY = $(addprefix media/,$(FIGS)) text.tex title.txt preamble.tex abstract.tex

all: manuscript.zip

manuscript.zip: $(TOCOPY) preprint.bbl preprint.tex $(STANDALONE)
	zip $@ $^

$(TOCOPY): %: ../% | media
	cp $^ $@

standalone-figures/fig-1.pdf: ../media/range-curves.pdf
standalone-figures/fig-2.pdf: ../media/s66-dists.pdf
standalone-figures/fig-3.pdf: ../_build/export-figures-figure0.pdf
standalone-figures/fig-4.pdf: ../_build/export-figures-figure1.pdf
standalone-figures/fig-5.pdf: ../media/scan-interp.pdf
preprint.bbl: ../_build/preprint.bbl
$(STANDALONE) preprint.bbl: | standalone-figures
	cp $^ $@

preprint.tex: ../preprint.tex
	sed 's;\\bibliography{refs-zotero,refs};\\input{preprint.bbl};' $^ >$@

standalone-figures media:
	mkdir $@

clean:
	-rm -r standalone-figures media title.txt preprint.bbl $(wildcard *.tex *.zip)
