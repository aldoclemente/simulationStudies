---
title: "Una breve introduzione ad RMarkdown"
subtitle: "A prova di Studente"
output: 
  pdf_document:
    number_sections: false
    fig_caption: true
    highlight: tango

urlcolor: blue
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.width= 4.5, 
                                   fig.height= 3.5, 
                                   fig.align= "center")
```

# Installazioni preliminari
1. Installare *RMarkdown*
```{r, eval=FALSE}
install.packages('rmarkdown')
```

2. Installare *TinyTex*
```{r, eval=FALSE}
install.packages("tinytex")
tinytex::install_tinytex()
```

# Crea il tuo primo file 
