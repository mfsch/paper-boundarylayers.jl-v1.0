# --- PATHS & FILES ---------------------------- #
tmpdir := .cache/tex
tooldir := .cache
pandocversion = 2.14
pandocbin = $(tooldir)/pandoc-$(pandocversion)/bin/pandoc
pandoccmd = $(pandocbin) "$<" --standalone --wrap=preserve --natbib
latexcmd = env TEXINPUTS=$(shell kpsewhich -var-value TEXINPUTS):../../../templates pdflatex --halt-on-error --interaction=batchmode
bibtexcmd = env BSTINPUTS="../../../templates:$(shell kpsewhich -var-value BSTINPUTS)" bibtex
lastrev := deb2d3e # for latex diffs


# --- PHONY TARGETS ---------------------------- #

.DEFAULT_GOAL := help
.PHONY: help default watch manuscript draft gmd diff
default: gmd
manuscript: manuscript-latest.pdf
gmd: manuscript-gmd.pdf
diff: manuscript-diff.pdf

help:
	@echo 'Run `make <goal>` with one of these goals:'
	@echo '- `manuscript`: Build manuscript PDF'
	@echo '- `gmd`: Build manuscript PDF with GMD template'
	@echo '- `watch`: Launch simulations to gather performance data'

watch:
	while inotifywait -e modify -r .; do make draft; done


# --- REFERENCES ------------------------------- #

biblio_path  := /home/mfsch/Notes/references
biblio_files := bibliography.bib
biblio_files_available := $(notdir $(wildcard $(addprefix $(biblio_path)/,$(biblio_files))))
ifneq ($(biblio_files),$(biblio_files_available))
    $(info -----------------------------------------------)
    $(info WARNING: Bibliography files missing!           )
    $(info Will not rebuild project-specific bibliography.)
    $(info -----------------------------------------------)
endif
bibtool_cmd = bibtool -q -s --preserve.key.case=on --delete.field={annote} --delete.field={date_added} --delete.field={file} --delete.field={isbn} --delete.field={note} --delete.field={number} --delete.field={url} --delete.field={tags}

references.bib: $(tmpdir)/manuscript.aux
	$(bibtool_cmd) -x $< $(addprefix -i $(biblio_path)/,$(biblio_files)) -o $@

$(tmpdir)/manuscript.bbl: $(tmpdir)/manuscript.aux references.bib
	cd $(<D) && $(bibtexcmd) $(basename $(notdir $<))

# --- DOCUMENTS -------------------------------- #

figures := $(wildcard ../figures/*.pdf)
$(tmpdir)/fig%: ../figures/fig%
	mkdir -p $(@D)
	cp $< $@

$(tmpdir)/manuscript.tex: manuscript.md ../templates/manuscript-draft.tex ../templates/outline.lua \
    $(patsubst ../figures/%,$(tmpdir)/%,$(figures)) | $(pandocbin)
	mkdir -p $(@D)
	$(pandoccmd) --output="$@" --metadata=bibliography:../../references.bib \
		--lua-filter="../templates/outline.lua" --lua-filter="../templates/linenomath.lua" \
		--template="../templates/manuscript-draft.tex" --metadata=graphics \
		"--metadata=date:$$(date '+%-d %b %Y')"

$(tmpdir)/manuscript.aux: $(tmpdir)/manuscript.tex
	cd $(<D) && $(latexcmd) --draftmode $(notdir $<)

manuscript-latest.pdf: $(tmpdir)/manuscript.tex $(tmpdir)/manuscript.bbl
	mkdir -p $(@D)
	cd $(dir $<) && $(latexcmd) --draftmode $(notdir $<) && $(latexcmd) $(notdir $<)
	mv "$(<:.tex=.pdf)" $@

draft: $(tmpdir)/manuscript.tex # only run one pass
	cd $(<D) && $(latexcmd) $(notdir $<)
	mv "$(<:.tex=.pdf)" manuscript-latest.pdf


# --- LATEX DIFF ------------------------------- #

$(tmpdir)/manuscript-$(strip $(lastrev)).tex:
	mkdir -p $(@D)
	git show $(strip $(lastrev)):./manuscript.md | \
		$(pandocbin) --output="$@" --standalone --wrap=preserve \
		--natbib --metadata=bibliography:../../references.bib \
		--lua-filter="../templates/outline.lua" \
		--template="../templates/manuscript-gmd.tex"

$(tmpdir)/manuscript-diff.tex: $(tmpdir)/manuscript-gmd.tex $(tmpdir)/manuscript-$(strip $(lastrev)).tex
	latexdiff --append-safecmd="mkern,url" --append-textcmd="codeavailability,codedataavailability,mkern" $(word 2,$^) $< > $@

manuscript-diff.pdf: $(tmpdir)/manuscript-diff.tex
	mkdir -p $(@D)
	cd $(<D) && $(latexcmd) --draftmode $(notdir $<)
	cd $(<D) && $(bibtexcmd) $(basename $(notdir $<))
	cd $(dir $<) && $(latexcmd) --draftmode $(notdir $<) && $(latexcmd) $(notdir $<)
	mv "$(<:.tex=.pdf)" $@


# --- LATEX JOURNAL TEMPLATES ------------------ #

$(tmpdir)/manuscript-%.tex: manuscript.md ../templates/manuscript-%.tex $(wildcard ../templates/*.lua) \
    $(patsubst ../figures/%,$(tmpdir)/%,$(figures)) | $(pandocbin)
	mkdir -p $(@D)
	$(pandoccmd) --output="$@" --metadata=bibliography:../../references.bib \
	    --lua-filter="../templates/outline.lua" --template="../templates/manuscript-$*.tex"
	cp $@ $@.bkp

manuscript-%.pdf: $(tmpdir)/manuscript-%.tex
	mkdir -p $(@D)
	cd $(<D) && $(latexcmd) --draftmode $(notdir $<)
	cd $(<D) && $(bibtexcmd) $(basename $(notdir $<))
	cd $(dir $<) && $(latexcmd) --draftmode $(notdir $<) && $(latexcmd) $(notdir $<)
	mv "$(<:.tex=.pdf)" $@


# --- SOFTWARE TOOLS --------------------------- #

$(pandocbin):
	@mkdir -p $(tooldir)
	wget -O $(tooldir)/pandoc.tar.gz https://github.com/jgm/pandoc/releases/download/$(pandocversion)/pandoc-$(pandocversion)-linux-amd64.tar.gz
	tar -xzf $(tooldir)/pandoc.tar.gz -C $(tooldir)
