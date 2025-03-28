---
title: 'PyKEP 3: A coolbox for interplanetary trajectory design'
tags:
  - Python
  - astrodynamics
  - interplanetary trajectories
  - preliminary design
  - optimization
authors:
  - name: Dario Izzo
    orcid: 0000-0002-9846-8423
    equal-contrib: false
    affiliation: "1," # (Multiple affiliations must be quoted)
  - name: Francesco Biscani
    equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: European Space Agency's Advanced Concepts Team, The Netherlands
   index: 1
 - name: European Space Agency's Advanced Concepts Team, The Netherlands
   index: 2
date: 28 March 2025
bibliography: paper.bib
---

# Summary

`PyKEP 3` is a Python toolbox developed at the [European Space Agency](https://www.esa.int) by the 
[Advanced Concpets Team](https://www.esa.int/act) to perform
quick analysis of interplanetary trajectory design problems. It is designed to be used by researchers
and engineers to prototype and test new ideas in the field of astrodynamics. The library provides
efficient implementations of algorithms for solving the multiple revolutions Lambert's problem, low-thrust
problems, multiple asteroid rendezvous problems, and more. It also provides support for [JPL SPICE](https://naif.jpl.nasa.gov/naif/toolkit.html),
SGP4 propagation, and the [Heyoka](https://bluescarni.github.io/heyoka.py/index.html) Taylor integration suite.

# Statement of need

`PyKEP 3` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from .... and support from ...

# References

