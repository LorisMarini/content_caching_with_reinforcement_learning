# Content caching with reinforcement learning
The problem of optimal data allocation in a network of wireless mobile terminals is known to be untractable for even small number of files and terminals (NP-Hard). This repository contains the code of the work published in IEEE Xplore: Distributed Caching based on Decentralized Learning Automata. 

# Problem 
Simply stated, the file placement problem or 'caching problem' arises when we want to find the optimal placement of F objects in H positions with a maximum of C objects per posiiton. Optimal refers to the allocation that minimizes some sort of cost function which in this context is the latency in the network. Trying all possible combination and permutation of the objects ('brute force' or 'exhaustive search' approach) becomes quickly unfeasible for small numbers of objects.

# Solution
There are many ways to approach a sub-optimal solution to the caching problem. We propose one that is inspired by a game of independent players (Learning Automata) which take actions and sense each others' choices to understand whether their strategy was good or not. Since there is no need of centralized entity that scores the player's choices, this approach is highly scalable. Under a simulated noisy environment our algorihtm approaches the performanc of a greedy strategy where every player minimizes their cost function. We propose a discrete generalized pursuit algorithm (DGPA

# Contributions

1) A modified version of discrete generalized pursuit algorithm (DGPA) based on the concept of conditional inaction, CI-DGPA with faster convergence speed and good accuracy.
2) A reward function that allows CI-DGPA to approach the performance of our benchmark (greedy algorithm)
3) Further convergence optimizations based on partitioning the search space
4) A draft of the comunication protocol needed in a realistic implementation of CI-DGPA in a wireless network.

# CODE

First run the following command: addpath(genpath('your-path/content_caching_with_reinforcement_learning')); where you have to be careful to replace 'your-path' with the folder path in your computer. The purpose of this is to tell matlab to add all directories and sub-directories to the path and know where to look for function definitions.

# 1. Network Deployment

To get started with the 'Configure_The_Network' function simply copy the text under the line "- EXAMPLE OF CALL -" in the function description and run it in the command window:

% Number of radio cells to deploy
Nrc = 2;
% Cell diameter [m]
Diam = 1000;
% Clear space around each helper [m]:
CS = 10;
% Step size
SS = 0.1
% Number of helpers in the network
Nh = 20;
% Number of users in the network
Nu = 100;
% A propagation Loss of 1 [km^-1]
Ap = 1;

[ Parameters ] = Configure_The_Network( Nrc, Diam, CS, SS, Nh, Nu, Ap)

You should e able to execute it without problems as well as debug the code. 


