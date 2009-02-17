(TeX-add-style-hook "main"
 (lambda ()
    (LaTeX-add-labels
     "fig:#3"
     "fig:#2")
    (TeX-add-symbols
     '("billedeorg" 3)
     '("billede" 4)
     "lodret")
    (TeX-run-style-hooks
     "hyperref"
     "colorlinks"
     "color"
     "usenames"
     "prettyref"
     "url"
     "calc"
     "textcomp"
     "amssymb"
     "amsmath"
     "shortvrb"
     "natbib"
     "square"
     "graphicx"
     "fontenc"
     "T1"
     "inputenc"
     "ansinew"
     "latex2e"
     "rep11"
     "report"
     "11pt"
     "a4paper"
     "forward"
     "abstract"
     "introduction"
     "theory"
     "evaluation"
     "conclusion"
     "biblio"
     "matlab_code")))

