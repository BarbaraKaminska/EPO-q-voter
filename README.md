# alpha-EPO-q-voter-model

This repository contains the code used to generate results for the publication "Thinking fast, acting slow: Self-anticonformity in a generalized expressed-private q-voter model". 

## The model 
State of each agent in the model is described by two binary variables: 
- S - expressed opinion (visible for other agents)
- &sigma; - private opinion (known only to the given voter)
Opinions &plusmn; 1 correspond to yes/no answers, being for/against given issue.
In the paper we also denote them using arrows, namely ↑ stands for +1 and ↓ for -1. 

Both opinions change in time as a result of independence or conformity. 
Figures below show algorithms of elementary changes. 

<img src="algorithm_harmony.png" width="500">
<img src="algorithm_dissonance.png" width="500">

### Model parameters
- N - number of nodes (agents in the system)
- q - size of influence group
- p - probability of independence
- &alpha; - probability of updating private opinion

## Monte Carlo simulations 
The C++ code allows the model to be run on a complete graph or an arbitrary network. 
To run the model on a complex network, the network must be provided in a .txt file as a list of neighbors. Each row is a list of neighbors of the following agents, which should be indexed from 0.

## Pair approximation 
The matlab codes allows to numerically solve pair approximation equations via evolution of system of ordinary differential equations until the stady state is reached. 
Note the notation u correspons to ↑ (up arrow), while d for ↓ (downarrow)
