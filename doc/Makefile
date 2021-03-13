all:		spline.pdf

spline.pdf:	spline.tex interpolate_avg.pdf positive_criteria.pdf spline.bbl spline.blg
		pdflatex spline.tex

interpolate_avg.pdf:	interpolate_avg.gp
		gnuplot interpolate_avg.gp

positive_criteria.pdf:	positive_criteria.gp
		gnuplot positive_criteria.gp

spline.bbl:	literature.bib spline.tex
		pdflatex spline; bibtex spline; pdflatex spline

spline.blg:	literature.bib spline.tex
		pdflatex spline; bibtex spline; pdflatex spline

clean:
		rm -f *.pdf *.aux *.dvi *.log *.ps *.bbl *.blg
