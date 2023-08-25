# Classification of Fano simplices

A julia package for the classification of and a database for Fano simplices. For background information see [arXiv:2308.12719](https://arxiv.org/abs/2308.12719).

## Description
The folders `fanosimp-d=2/` to `fanosimp-d=4/` contain text file, one per Gorenstein index.
For instance the file `fanosimp-d=2-g=5.txt` contains the Fano simplices of dimension two and Gorenstein index five.
Each line of one of these files is a simplex, presented as a list of vertices. For instance the entry `[[1, 0], [1, 2], [-6, -5]]`
is the Fano triangle with vertices $(1,0)$, $(1,2)$ and $(-6,-5)$.
