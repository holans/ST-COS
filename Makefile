PROJNAME := stcos
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_path := $(dir $(mkfile_path))

## See http://kbroman.org/pkg_primer/pages/docs.html
## for help on generating Roxygen2 documentation
## NOTE: Building the pdf seems not to work until the package is already installed

pdf: manual
	R CMD Rd2pdf --force --no-preview --output=$(PROJNAME).pdf .

manual:
	R -e 'library(devtools); document(roclets = c("collate", "rd"))'

clean:
	rm -rf man $(PROJNAME).pdf src/$(PROJNAME).so
	find $(current_path)/src -type f -name '*.o' -exec rm {} \;
