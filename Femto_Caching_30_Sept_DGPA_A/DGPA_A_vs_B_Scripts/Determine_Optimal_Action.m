function [ This_Optimum ] = Determine_Optimal_Action( Real_Optimum, Variance,S )

%   'Real_Optimum' is the position of the optimal action in the space S of
%   possible actions.'Variance' represents the noise of what the
%   environment considers to be optimal at a given time. 
%   Typical Values for Variance are 1,2,3,5,10;

% Define the array of actions
V = 1:1:length(S);

% Extract the probability density function for a Gaussian with 
Prob = pdf('Normal', V, Real_Optimum, Variance);

% Normalize the probability
Norm_Prob = Prob./sum(Prob);

% Use the normalized probability density function 

This_Optimum = randsrc(1,1,[V;Norm_Prob]);

end
