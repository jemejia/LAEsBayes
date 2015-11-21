NAME  = laes_bayes
LATEX = pdflatex
VIEWER = evince
TEXEDIT = emacs

$(NAME).pdf:	$(NAME).tex  makefile references.bib
		$(LATEX) $(NAME).tex 
		bibtex ${NAME}
		$(LATEX) $(NAME).tex 
		$(LATEX) $(NAME).tex 
clean:
		rm -f $(NAME).aux
		rm -f $(NAME).out
		rm -f $(NAME).bbl
		rm -f $(NAME).log

view: 
		$(VIEWER) $(NAME).pdf &

