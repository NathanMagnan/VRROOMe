# VRROOMe

Vector Resonant Relaxation Out Of Maximum entropy

## Installation

You'll need to install Julia, by following the instructions at ` https://julialang.org/downloads/platform/ `\
To invoke Julia in the terminal, you'll need to make sure that the Julia command-line program is in your ` PATH `

You'll also need to install the following Julia packages:
- Distributions
- Plots
- Interpolations
- HDF5
- StaticArrays
- FastGaussQuadrature
- Random
- LinearAlgebra

**WARNING: Do NOT interupt the installation of the packages!**

Once all of this is done, you can test if everything is ready by running the file ` Code/Example.jl ` in the REPL.\
After a minute, it should produce one line of warning, don't worry about it:

> ` WARNING: using Distributions.scale in module Main conflicts with an existing identifier. `

It should also display 2 plots, which should be perfect replicas of ` Figures/Initial_cluster.png ` and ` Figures/Relaxed_cluster.png `.

## First steps

The ` Code/Example.jl ` has been designed to provide examples of the commands I find the most usefull.\
Please read through it to get started.

## Advanced usage

The more advanced reader might feel the need to play with lower-level functions. Please do so, they are all in the ` Code/ ` folder.\
I tried to document these within the code itself. But should my comments be unsufficient, please send me an e-mail at ` nathan.magnan@maths.cam.ac.uk ` to request clarifications.