function [ LA, I ] = CI_DGPA( Actions_Space ,Optimal_Action ,Environment_Variance ,Resolution ,Initial_Iterations )
% DGPA_Signle_Player 
% This function performs a search of the optimal action in a space of actions S, 
% based on the Discrete Generalised Pursuit Algorithm for Learning
% Automata.

% THis is version B. With a Different Approach in the Probability updating
% policy. Action vector a(t) defined as follows:

% FOR j ~= i
% aj(t) = D, if dj(t) =< di(t) && Pj ~= 0
% aj(t) = I, if dj(t) > di(t) 
% aj(t) = N, if Pj == 0

% FOR j = i
% aj(t) = I, if di(t) == max(dj(t))
% aj(t) = D, otherwise.

% D = Decrement;
% I = Increment;
% N = Don't do anything.

% INPUTS
% 'S' should be a vector of integers from 1 to length(S) so that the 5th element of S is action  number 5.
% 'Optiomal_Action' = Expectation of the Gaussian Random Environment; 
% 'Environment_Variance' =  Variance of the Envoronment (typical values are 0.5,1,2,3).
% 'Initial_Iterations' should be a small number, generally comprised betwee 10 and 100.
% 'N' = Resolution. Determines the width of the discrete increments and decrements of the probability vector P at each iteration. 
% (Higher 'N' leads to longer convergence times but higher accuracy). 

% OUTPUTS
% The function returns the Learned Action 'LA' as well as the number of
% Iterations 'I' performed to learn.

% Author: Loris Marini 
% Version: 1.0.0 
% Date: 03/09/2014

%% VARIABLES INITIALISATION

    S = Actions_Space;  % Space of Actions
    r = length(S);    % Number of actions to search from
    P = (1/r).*ones(1,r);  % Actions Probability Vector uniformally initialised
    D = zeros(1,r);   % Reward Probability Estimates
    W = zeros(1,r);   % Incidence of Rewards 
    Z = zeros(1,r);   % Actions Selection Incidence
    N = Resolution;   % Resolution Parameter
    Delta = 1/(r.*N); % Resulution Step
    N_Ini = Initial_Iterations; % Number of routines for initialisation of estimates D
    Opt_Act = Optimal_Action;      % The REAL OPTIMAL ACTION
    Env_Variance = Environment_Variance; 

    %% Realisation of the Normal Random Environment
    % We first determine what is the optimal action for the environment now.
    % Current_Optimum = Determine_Optimal_Action(Opt_Act, r/10, S);

    
     %% MLE Initialisation of the whole vector D

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

    % [D,W,Z] = MRC_Initialisation( Z,W,S,P,Opt_Act,Env_Variance,N_Ini);

    %% Variables allocation and Initialisation

    Pre_All = 1000;

    % Variables Allocation
    New_P = zeros(Pre_All, r);
    New_D = zeros(Pre_All, r);
    K = zeros(1,Pre_All);
    Increment = zeros(1,Pre_All);
    Decrement = zeros(1,Pre_All);

    Post_Inc_Rem_Sum = zeros(1,Pre_All);
    Pre_Inc_Rem_Sum = zeros(1,Pre_All);

    Excess = zeros(1,Pre_All);
    Consistent_Increment = zeros(1,Pre_All);
    New_Delta = zeros(1,Pre_All);

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

     i = 2;

    while ~sum ( New_P(i-1,:) > 0.999 )

        Current_Optimum = Determine_Optimal_Action(Opt_Act, Env_Variance, S);

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
        
        % 3) Update the reward probability estimates
        New_D(i,:) =  W ./ Z; 
        
        % 4) Determine which action should be incremented and which decremented.
        Direction_Vector = zeros(1,r);
        % Vector element = 1  INCREMENT
        % Vector Element = 0  DECREMENT
        % Vector Element = Inf Don't do anything.
        
        Actions_Subset = New_P( i-1,:) ~= 0;
        
        for o = 1:1:r
            if (o ~= Ai)
                if ( New_D(i-1,o) <= New_D(i-1,Ai)  &&  New_P( i-1, o )~= 0) 
                    Direction_Vector(o) = 0;
                elseif ( New_D(i-1,o) > New_D(i-1,Ai) &&  New_P( i-1, o )~= 0)                      % ( New_D(i-1,o) > New_D(i-1,Ai)  &&  New_P( i-1, o )~= 0)
                    Direction_Vector(o) = 1;
                elseif ( New_P( i-1, o ) == 0)
                    
                    Direction_Vector(o) = Inf;
                end      
            elseif (o == Ai)
                [Value, Position] = max(New_D(i-1,Actions_Subset));
                Number_of_Maxima = sum ( New_D(i-1,Actions_Subset) == Value );
                 if Number_of_Maxima == 1    % Single Maximum
                    if (New_D(i-1,Ai) == Value)
                        Direction_Vector(o) = 1;
                    elseif (New_D(i-1,Ai) ~= Value)
                        Direction_Vector(o) = 0;
                    end
                 elseif Number_of_Maxima > 1  % Multiple Maxima
                     Maxima =  New_D(i-1,:) == Value;
                     Index = Maxima & Actions_Subset;
                     if sum ( New_D(i-1,Index) == New_D(i-1,Ai) ) > 1
                          Direction_Vector(o) = 1;
                     else
                          Direction_Vector(o) = 0;
                     end
                 end
            end
        end
                
        K(i) = sum( Direction_Vector == 1 );
        N_Dec(i) = sum( Direction_Vector == 0 );
        Inactions(i) = sum(Direction_Vector == Inf);
        if (K(i) + N_Dec(i) +  Inactions(i)) ~= r
            error('Inconsistency.');
        end
        Who_To_Inc =  Direction_Vector == 1;
        Who_To_Dec =  Direction_Vector == 0;
        
        if (K(i) == 0)
            error('There should always be an action to make more likely...');
        end
        
        % 2. Determine consistent Increments and Decrements.
        if ( K(i)~=0 )
            Increment(i) = Delta./ K(i);
        else
            Increment(i) = 0;
        end
        if ( N_Dec(i) ~=0)
            Decrement(i) = - ( Delta ./ N_Dec(i) );
        else
            Decrement(i) = 0;
        end
        
    %% NEW P FOR ACTIONS MORE LIKELY TO BE REWARDED   
              
        New_P( i, Who_To_Inc ) = New_P( i-1, Who_To_Inc ) + Increment(i);
        
        % CHECK CONVERGENCE 1
        learned = S( New_P(i,:) > 0.999);
        if (learned)
            I=i;
            LA = learned;
            if ( learned == Opt_Act)
                disp (['CC1: The algorithm converged to the CORRECT action number ' num2str( S(learned)) ] );
                disp (['CC1: Total number of iterations ' num2str(i)] );
                break
            else
                disp (['CC1: The algorithm converged to the w.r.o.n.g. action number ' num2str( S(learned)) ] );
                disp (['CC1: Total number of iterations ' num2str(i)] );
                break
            end
        end
        
        New_P( i, Who_To_Dec ) = max ( New_P( i-1, Who_To_Dec ) + Decrement(i), zeros(1,N_Dec(i)) );
        
        if (sum( New_P( i, : )) > 1)
            
            Excess(i) = abs(sum( New_P( i, : )) - 1 );
            Penalty_Decrements(i) = Excess(i)/ K(i);
            New_P( i, Who_To_Inc ) = New_P( i, Who_To_Inc ) -  Penalty_Decrements(i);  
        end
        
         if (sum( New_P( i, : )) > 1+10^-12 )
            error('It didnt work. ERROR 101');
         end
         if (sum( New_P( i, : )) < 1-10^-12 )
             error('It didnt work. ERROR 102');
         end
          
       %% LAST CONSISTENCY CHECK
        % 7) 
        who_is_negative = New_P(i,:) < 0;
        N_negatives = sum(who_is_negative);
        Sum_Probability = sum (New_P(i,:));

        if ( N_negatives ~= 0 ) 
            error('The error occured at iteration number %d. We had %d negative probability(ies). The remaining sum was %d',i,N_negatives,Remaining_Sum(i));
        end    
        if (Sum_Probability > (1 + 10^-12) )
           error('The actions probabilities at iteration %d add up to %d which is greater than 1.', i,Sum_Probability);
        end

        i = i+1;
      
    end  % This closes the while loop.
    
 
end
