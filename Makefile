objects := $(wildcard R/*.R) DESCRIPTION $(wildcard src/*.cpp)
version := $(shell egrep "^Version" DESCRIPTION | awk '{print $$NF}')
pkg := $(shell egrep "^Package" DESCRIPTION | awk '{print $$NF}')
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log
no_tex_flag := --no-build-vignettes --no-manual


.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

$(tar): $(objects)
	rm -rf R/RcppExports.R src/RccpExports.cpp
	Rscript -e "library(methods); Rcpp::compileAttributes();" \
	-e "devtools::document();";
	R CMD build $(no_tex_flag) .

$(checkLog): $(tar)
	R CMD check $(no_tex_flag) $<

## make tags
.PHONY: TAGS
TAGS:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"

## render readme file
.PHONY: readme
readme: README.md
README.md: README.Rmd
	@Rscript -e "rmarkdown::render('$<')"

## do some cleaning
.PHONY: clean
clean:
	@rm -rf *~ */*~ *.Rhistroy src/{*.o,*.so} *.tar.gz *.Rcheck/ .\#*
