#!/bin/sh

name=poptest

if [ x"$1" = x"clean" ]; then
	rm -f $name.aux $name.dvi $name.log $name.toc \
		$name.nav $name.snm $name.toc \
		$name.bbl $name.blg $name.out
	exit 0
fi

#pdflatex $name.tex || exit 1
#bibtex $name.aux || exit 1
pdflatex $name.tex && pkill -HUP mupdf
