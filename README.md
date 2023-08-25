# Classification of Fano simplices

A julia package for the classification of and a database for Fano simplices. For background information see [arXiv:2308.12719](https://arxiv.org/abs/2308.12719).

## Description
The folders `fanosimp-d=2/` to `fanosimp-d=4/` contain text file, one per Gorenstein index.
For instance the file `fanosimp-d=2-g=5.txt` contains the Fano simplices of dimension two and Gorenstein index five.
Each line of one of these files is a simplex, presented as a list of vertices. For instance the entry `[[1, 0], [1, 2], [-6, -5]]`
is the Fano triangle with vertices $(1,0)$, $(1,2)$ and $(-6,-5)$.

## How to use `fanosimp.jl`
The data was produced with the julia package `fanosimp.jl`. There is no separate installation, just download the file and include it in your script.
fanosimp uses the following packages:
- `Combinatorics`
- `Primes`
- `AbstractAlgebra`

These have to be installed before using fanosimp.

Navigate to the folder that contains the file `fanosimp.jl` and start a Julia instance. The classification routine is multithreaded. To make use of this, start Julia with
```
> julia -t auto
```
to use all available threads, or type
```
> julia --threads=4
```
to restrict to a custom number of threads, in this case four.

Let's say we want to classify the reflexive triangles. We start julia with four threads, include the fanosimp package and initialize the classification. Here is what the output looks like:
```
> julia --threads=4
julia> include("fanosimp.jl");
julia> fullclass_fanosimp(2,1)
Classification of Fano simplicies of dimension 2 and Gorenstein index 1.
Using 4 threads.
Computing list of weight systems...
Done. List of 3 well-formed weight systems.
Starting classification. This might take a while.
[1, 1, 1] computing...
[2, 1, 1] computing...
[3, 2, 1] computing...
[3, 2, 1] done. 1/3
[2, 1, 1] done. 2/3
[1, 1, 1] done. 3/3
Done.
Found 5 Fano simplices of dimension 2 and Gorenstein index 1.
```
**The procedure fullclass_fanosimp creates a folder `temp/` in the current working directory**. vIn this folder there are now three text files, corresponding to the three weight systems, containing the reflexive triangles.
