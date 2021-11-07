# VRROOMe

Vector Resonant Relaxation Out Of Maximum entropy

## Installation

You'll need to install Julia, by following the instructions at ` https://julialang.org/downloads/platform/ `\
To invoke Julia in the terminal, you'll need to make sure that the Julia command-line program is in your ` PATH `

You'll also need to install the following Julia packages:
- Distributions
- Plots
- PlotlyJS
- Interpolations
- HDF5
- StaticArrays
- FastGaussQuadrature
- Random
- LinearAlgebra

**WARNING: Do NOT interupt the installation of the packages!**

## Testing the installation with the example script

Once all of this is done, you can test if everything is ready by running the file ` Code/Example.jl ` in the REPL.\
After a few minutes, it should produce:

> ` Starting -- Julia seems to be working. `\
> \
> ` Let's first create an initial cluster, with 16 star formation events `\
> ` WARNING: using Distributions.scale in module Main conflicts with an existing identifier. `\
> ` We can compute this cluster's binding energy E and spin s `\
> ` The cluster has been created and its (E, s) computed `\
> ` Now let's plot it `\
> \
> ` Now let's create the initial DF related to this initial cluster `\
> ` WARNING: redefinition of constant TAB_INT_SL. This may fail, cause incorrect answers, or produce other errors. `\
> ` The initial DF related to the initial cluster has been created `\
> \
> ` Let's modify this DF so that it maximizes the entropy S `\
> ` true `\
> ` The entropy has been maximized `\
> ` Finally, let's plot the final DF `

Don't worry about the warnings.\
It should also display 2 plots, which should be perfect replicas of ` Figures/Initial_cluster.png ` and ` Figures/Relaxed_cluster.png `:

` Figures/Initial_cluster.png ` | ` Figures/Relaxed_cluster.png `
-|- 
<img src="Figures/Initial_cluster.png" alt="drawing" width="400"/> | <img src="Figures/Relaxed_cluster.png" alt="drawing" width="400"/>

If you prefer to run Julia from the terminal rather than the REPL, you might get 2 blank windows instead of 2 plots.\
This is a registered bug in PlotlyJS v0.14.1, you have to upgrade to version v0.18.8 or higher:

> ` Pkg.update("PlotlyJS") `

Please note that from terminal, the plots might take some time to appear, even after the window pops up.

## First steps

The ` Code/Example.jl ` has been designed to provide examples of the commands I find the most usefull.\
Please read through it to get started.

## Advanced usage

The more advanced reader might feel the need to play with lower-level functions. Please do so, they are all in the ` Code/ ` folder.\
I tried to document these within the code itself. But should my comments be unsufficient, please send me an e-mail at ` nathan.magnan@maths.cam.ac.uk ` to request clarifications.

## Usage in research

Should you use VROOMe for research, we invite you to cite the paper **To Be Completed once it is published**.
