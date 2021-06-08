.PHONY: all
rmd := $(patsubst %.R,%.Rmd,$(wildcard *.R)) \
	   $(patsubst %.py,%.Rmd,$(wildcard *.py))  \
	   $(patsubst %.sh,%.Rmd,$(wildcard *.sh))

.PRECIOUS: %.Rmd

outmd := $(patsubst %.Rmd,%.out.md,$(rmd))
outhtml := $(patsubst %.Rmd,%.html,$(rmd))

html: $(outhtml)
all: $(rmd)
md: $(outmd)

.PHONY: clean
clean:
	@echo rm -rf *_cache/ *_files/ *.html

%.Rmd: %.R
	jupytext --to Rmd '$^'

%.Rmd: %.py
	jupytext --to Rmd '$^'

%.Rmd: %.sh
	jupytext --to Rmd '$^'

%.md %.html: %.Rmd
	Rscript -e 'rmarkdown::render("$^", rmarkdown::html_document(keep_md=T), output_file="$@")'
