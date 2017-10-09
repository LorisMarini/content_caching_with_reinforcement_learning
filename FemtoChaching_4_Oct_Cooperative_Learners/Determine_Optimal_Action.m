function [ This_Optimum ] = Determine_Optimal_Action( Real_Optimum, Variance,S )
%UNTITLED Summary of this function goes here
%   'Real_Optimum' is the position of the optimal action in the space S of
%   possible actions.'Variance' represents the variance of what the
%   environment considers to be optimal at a given time. 

%   Typical Values for Variance are 1,2,3,5,10;

V = 1:1:length(S);
Prob = pdf('Normal', V, Real_Optimum,Variance);
Norm_Prob = Prob./sum(Prob);

This_Optimum = randsrc(1,1,[V;Norm_Prob]);

end
