% This script performs a search of the optimal action in a space of actions S, 
% based on the Discrete Generalised Pursuit Algorithm for Learning
% Automata.

% Author: Loris Marini 
% Version: 1.0.0 
% Date: 20/08/2014

%% VARIABLES INITIALISATION
clear all
close all
S = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];  % Space of Actions
r = length(S); % Number of actions to search from
P = (1/r).*ones(1,r);  % Actions Probability Vector uniformally initialised
D = zeros(1,r);   % Reward Probability Estimates
W = zeros(1,r);   % Rewards Incidence 
Z = zeros(1,r);   % Actions Selection Incidence
N = 1;            % Resolution Parameter
Delta = 1/(r.*N); % Resulution Step
N_Ini = 100;      % Number of routines for initialisation of estimates D
Opt_Act = 9;      % The REAL OPTIMAL ACTION
N_Iterations = 10;

%% Realisation of the Normal Random Environment
% We first determine what is the optimal action for the environment now.
% Current_Optimum = Determine_Optimal_Action(Opt_Act, r/10, S);

 %% MLE Initialisation of the whole vector D
 
 while  sum( Z < N_Ini ) ~= 0
     
     Ai = randsrc(1,1,[S;P]);
     Z(Ai) = Z(Ai) + 1;
     
     Current_Optimum = Determine_Optimal_Action(Opt_Act, r/10, S);
     
     if (Environment(Ai,S,Current_Optimum))       
         W(Ai) = W(Ai) + 1;
     else
         W(Ai) = W(Ai);
     end
 end

 D = W./Z;
     
%% Variables allocation and Initialisation

% Variables Allocation
New_P = zeros(N_Iterations, r);
New_D = zeros(N_Iterations, r);
K = zeros(1,N_Iterations);
Increment = zeros(1,N_Iterations);
Decrement = zeros(1,N_Iterations);
Remaining_Sum = zeros(1,N_Iterations);
Excess = zeros(1,N_Iterations);
Consistent_Increment = zeros(1,N_Iterations);
New_Delta = zeros(1,N_Iterations);

% Variables Initialisation 
New_P (1,:) = P;
New_D (1,:) = D;
Increment(1) = 0;
Decrement(2)= 0;
K(1) = 0;
Remaining_Sum(1) = 0;
New_Delta(:) = Delta;

%% DGPA at work to learn the optimum action.

 % 1) We first select an action according to the action probability vector
 % 2) We collect the Feedback from the Environment and update the reward vector
 % 3) We determine how many action should be made more likely
 % 4) We determine Increments and Decrements for the actions probabilities
 % 5) We Update the action's probability vector P
 % 6) We update the estimates of the reward probabilities
 % 7) Error Check.

for i=2:1:N_Iterations 
    
    Current_Optimum = Determine_Optimal_Action(Opt_Act, r/10, S);
    
    % 1)
    Ai = randsrc(1,1,[S; New_P(i-1,:)]);      % Random selection of action in position Ai
    Z(Ai) = Z(Ai) + 1;                        % Update the counter for action number Ai
    
    % 2)
    Feedback = Environment(Ai, S, Current_Optimum);
    if (Feedback)
        W(Ai) = W(Ai) + 1;
    else
        W(Ai) = W(Ai);
    end
    
   % 3) 
    K(i) = sum( New_D(i-1,:) >  New_D(i-1,Ai));
    
   % 4)   
   if (K(i)~=0)
       Increment(i) = Delta./K(i);
   else
       Increment(i) = 0;
   end
   if ((r - K(i)) ~=0)
       Decrement(i) = - (Delta ./ (r - K(i)));
   else
       Decrement(i) = 0;
   end
    
    Who_To_Inc = New_D(i-1,:) > New_D(i-1,Ai);
    Who_To_Dec = New_D(i-1,:) < New_D(i-1,Ai);
    N_Inc = sum( Who_To_Inc );
    N_Dec = sum( Who_To_Dec );
    
    % Will the elements who are going to be incremented add up to more than 1 ?
    % We need to guarantee consistency. Therefore, the maximum allowable increment 
    % should not cause the sum of incremented probabilities to exceed 1.
    % We should change Delta, because it appears also in the decrements.
    
    % 5)

    %%  NEW P FOR ACTIONS LESS LIKELY TO BE REWARDED
    
    if N_Dec ~= 0
        New_P( i, Who_To_Dec ) = max ( New_P( i-1, Who_To_Dec ) + Decrement(i), zeros(1,N_Dec) );
    end
    
     %% NEW P FOR ACTIONS MORE LIKELY TO BE REWARDED
     
    if N_Inc ~= 0
        New_P( i, Who_To_Inc ) = min ( New_P( i-1, Who_To_Inc ) + Increment(i), ones(1,N_Inc) );
        Remaining_Sum(i) =  sum(New_P( i, 1:Ai-1)) + sum( New_P( i, Ai+1:end) );
        
        if Remaining_Sum(i) > 1           
            % The proposed increment causes an excess. It has to be adjusted.
            Excess(i) = abs(1- Remaining_Sum(i) );
            New_Delta(i) = New_Delta(1) - Excess(i);
            Consistent_Increment(i) =  max( New_Delta(i)./K(i), 0 );
            New_P( i, Who_To_Inc ) = min ( New_P( i-1, Who_To_Inc ) + Consistent_Increment(i), ones(1,N_Inc) );
        end 
    end
    
    %%  NEW P FOR CURRENTLY CHOOSEN ACTION
    
    Remaining_Sum(i) =  sum(New_P( i, 1:Ai-1)) + sum( New_P( i, Ai+1:end) );
       
    if (Remaining_Sum(i) > 1)
        if (abs(1 -  Remaining_Sum(i)) < 10^-12)
            
            New_P(i, Ai) = 0;
        else
            error('What should I do now?');
        end
    else
        New_P(i, Ai) = 1 -  Remaining_Sum(i);
    end
    
    %% UPDATE THE REWARD ESTIMATES
    % 6)
    New_D(i,:) =  W ./ Z; 
   
    %% CONSISTENCY CHECK
    % 7) 
    who_is_negative = New_P(i,:) < 0;
    N_negatives = sum(who_is_negative);
    Sum_Probability = sum (New_P(i,:));
    
    if ( N_negatives ~= 0 ) 
        error('The error occured at iteration number %d. We had %d negative probability(ies). The remaining sum was %d',i,N_negatives,Remaining_Sum(i));
    end    
    if (Sum_Probability > (1 + 10^-12) )
       error('The actions probabilities at iteration %d add up to %d.', i,Sum_Probability);
    end
    
   
end

stop =1;
