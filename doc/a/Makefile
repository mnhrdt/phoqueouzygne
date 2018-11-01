all: a.pdf #p.pdf

%.pdf: %.tex
	grep ^%SCRIPT $< | cut -c9- | sh
	latexmk $<

clean:
	$(RM) o/*png o/*pdf
	latexmk -quiet -C a p
