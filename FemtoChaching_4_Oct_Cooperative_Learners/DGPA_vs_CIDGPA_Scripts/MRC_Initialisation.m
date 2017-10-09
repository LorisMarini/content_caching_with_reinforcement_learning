function [ D,W,Z ] = MRC_Initialisation( Z,W,S,P,Opt_Act,Env_Variance,N_Ini)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

while  sum( Z < N_Ini ) ~= 0
     
     Ai = randsrc(1,1,[S;P]);
     Z(Ai) = Z(Ai) + 1;
     
     Current_Optimum = Determine_Optimal_Action(Opt_Act, Env_Variance, S);
     
     if (Environment(Ai,S,Current_Optimum))       
         W(Ai) = W(Ai) + 1;
     else
         W(Ai) = W(Ai);
     end
end

 D = W./Z;
 
end

