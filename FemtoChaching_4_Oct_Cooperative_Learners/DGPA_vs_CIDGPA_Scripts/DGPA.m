function [ LA, I ] = DGPA( Actions_Space ,Optimal_Action ,Environment_Variance ,Resolution ,Initial_Iterations )
% DGPA_Signle_Player 
% This function performs a search of the optimal action in a space of actions S, 
% based on the Discrete Generalised Pursuit Algorithm for Learning
% Automata.

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
% Date: 21/08/2014

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
            Pre_Inc_Rem_Sum(i) =  sum(New_P( i, 1:Ai-1)) + sum( New_P( i, Ai+1:end) ); 
          
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
            
           %% INCREMENT CONSISTENCY CHECK         
            if Pre_Inc_Rem_Sum(i) > 1           
              
                Excess(i) = abs(1- Pre_Inc_Rem_Sum(i) );
                New_Delta(i) = New_Delta(1) - Excess(i);
                Consistent_Increment(i) =  max( New_Delta(i)./K(i), 0 );
                New_P( i, Who_To_Inc ) = min ( New_P( i-1, Who_To_Inc ) + Consistent_Increment(i), ones(1,N_Inc) );

                % CHECK CONVERGENCE 2
                learned = S( New_P(i,:) > 0.999); 
                if (learned) 
                    I=i;
                    LA = learned;
                    if ( learned == Opt_Act)
                        disp (['CC2: The algorithm converged to the CORRECT action number ' num2str( S(learned)) ] );
                        disp (['CC2: Total number of iterations ' num2str(i)] );
                        break
                    else
                        disp (['CC2: The algorithm converged to the w.r.o.n.g. action number ' num2str( S(learned)) ] );
                        disp (['CC2: Total number of iterations ' num2str(i)] );
                        break
                    end
                end
                
            end
        end

        %%  NEW P FOR CURRENTLY CHOOSEN ACTION
        
        Post_Inc_Rem_Sum(i) =  sum(New_P( i, 1:Ai-1)) + sum( New_P( i, Ai+1:end) );
        
        if ( Post_Inc_Rem_Sum(i) < 1 )
            
            New_P(i, Ai) = 1 -  Post_Inc_Rem_Sum(i);
                          
        elseif ( Post_Inc_Rem_Sum(i) > 1 )
            if (abs(1 -  Post_Inc_Rem_Sum(i)) < 10^-12)
                
                New_P(i, Ai) = 0;
            else
                error('There was a mistake in the probability computation. Check the code.');
            end
        end
        
        % CHECK CONVERGENCE 3
        learned = S( New_P(i,:) > 0.999); 
        if (learned)
            I=i;
            LA = learned;            
            if ( learned == Opt_Act)
                disp (['CC3: The algorithm converged to the CORRECT action number ' num2str( S(learned)) ] );
                disp (['CC3: Total number of iterations ' num2str(i)] );
                break
            else
                disp (['CC3: The algorithm converged to the w.r.o.n.g. action number ' num2str( S(learned)) ] );
                disp (['CC3: Total number of iterations ' num2str(i)] );
                break
            end
        end
                
            
        %% UPDATE THE REWARD ESTIMATES
        % 6)
        New_D(i,:) =  W ./ Z; 

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
