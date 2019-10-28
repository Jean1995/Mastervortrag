all: build/mastervortrag.tex


texoptions = \
	     --lualatex \
	     --interaction=nonstopmode \
	     --halt-on-error \
	     --output-directory=build

build/mastervortrag.tex: FORCE | build
	latexmk $(texoptions) mastervortrag.tex

preview: FORCE | build
	latexmk $(texoptions) -pvc mastervortrag.tex

FORCE:

build:
	mkdir -p build

clean:
	rm -r build
