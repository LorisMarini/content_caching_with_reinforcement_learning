function [ LA, I ] = DGPA( Actions_Space ,Optimal_Action ,Environment_Variance ,Resolution ,Init_Iterations )
 
%{
-------------------------   AUTHORSHIP  -------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

-------------------------   DESCRIPTION   -------------------------

This function performs a search of the optimal action in a space of actions S, 
based on the Discrete Generalised Pursuit Algorithm for Learning Automata.
This is a general function and does not care what the action space actually
represents. 

------------------------- INPUT PARAMETERS -------------------------

-- Actions_Space -- 
   A 1D array of integers from 1 to length(S) so that the i-th element of S is action number i.

-- Optimal_Action -- 
   Expectation of the Gaussian Random Environment. That is the expected
   value of what the environment considers optimal action.
       
-- Environment_Variance --
   Variance of the Envoronment (typical values are 0.5,1,2,3). It
   represents the noise in the system, how uncertain is the optimal action
   at any given time. 

-- Resolution --
   Determines the width of the discrete increments and decrements of the 
   probability vector P at each iteration. Is related to \Delta in
   eq.(10) and eq.(14) of Marini, L., Li, J., & Li, Y. (n.d.). "Distributed 
   Caching based on Decentralized Learning Automata, 1–6" via Delta = 1/(r.*N). 
   Higher Resolution leads to longer convergence times but higher accuracy.
    
-- Init_Iterations --
   Is the numebr of times (during initialization) that actions must be selected 
   before temrinating the initialization. Should be comprised betwee 10 and 100.


------------------------- OUTPUT PARAMETERS -------------------------

-- LA -- 
    The learned action from the set.
     
-- I --
    The number of iterations that were needeed to learn.

------------------------- EXAMPLE OF CALL -----------------------

Actions_Space = 1:1:15;
Optimal_Action = 7;
Environment_Variance = 3;
Resolution = 100;
Init_Iterations = 100;
[ LA, I ] = DGPA_Single_Player( Actions_Space ,Optimal_Action ,Environment_Variance ,Resolution ,Init_Iterations )

% ----------------------------   CODE     --------------------------
%}


%% Variable Initialization

    S = Actions_Space;                   % Space of Actions
    r = length(S);                       % Number of actions to search from
    P = (1/r).*ones(1,r);                % Actions Probability uniformally initialised
    D = zeros(1,r);                      % Reward Probability Estimates
    W = zeros(1,r);                      % Number of times actions are rewarded 
    Z = zeros(1,r);                      % Number of times actions are selected 
    N_Ini = Init_Iterations;             % Number of routines for initialisation of estimates D
    Opt_Act = Optimal_Action;            % The real optimal action
    Env_Variance = Environment_Variance; % the noise in the environment. 
    
    N = Resolution;                      % Resolution Parameter.
    Delta = 1/(r.*N);                    % Resulution Step. 

    % N is simply an integer. Delta is the \Delta in eq.(10) and eq.(14) of
    % the paper.
    
    %% Maximum Likelihood estimation (MLE) of the vector D
    
    [D,W,Z] = MRC_Initialisation( Z,W,S,P,Opt_Act,Env_Variance,N_Ini);
    

    %% Variables allocation and Initialisation

    Pre_All = 1000;

    % Variables pre-allocation
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

    %% DGPA at work to learn the optimal action.

     % 1) Select an action according to the action probability vector
     % 2) Collect the Feedback from the Environment and update the reward vector
     % 3) Determine how many action should be made more likely
     % 4) Determine amount of Increments and Decrements for the actions probabilities
     % 5) Update the action's probability vector P
     % 6) Update the estimates of the reward probabilities
     % 7) Error Check
     % 8) Repeat from 1 until a certai threshold is met
     % 9) Terminate the execution of the algorithm and assume learning complete.

     i = 2;

    while ~sum ( New_P(i-1,:) > 0.999 )

        Current_Optimum = Determine_Optimal_Action(Opt_Act, Env_Variance, S);

        % 1) Random selection of action Ai
        Ai = randsrc(1,1,[S; New_P(i-1,:)]);
        %  Update the counter for action Ai
        Z(Ai) = Z(Ai) + 1;                        

        % 2) collect the Feedback from the Environment and update the reward vector
        Feedback = Environment(Ai, S, Current_Optimum);
        if (Feedback)
            W(Ai) = W(Ai) + 1;
        else
            W(Ai) = W(Ai);
        end

       % 3) We determine how many action should be made more likely, that is the K(t)
       % in eq(10) & eq.(14) of Marini, L., Li, J., & Li, Y. (n.d.). "Distributed 
       % Caching based on Decentralized Learning Automata, 1–6".
       
        K(i) = sum( New_D(i-1,:) >  New_D(i-1,Ai));

       % 4) determine Increments and Decrements for the actions probabilities
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

        % 5) Update the action's probability vector P
        %    FOR ACTIONS LESS LIKELY TO BE REWARDED

        if N_Dec ~= 0
            New_P( i, Who_To_Dec ) = max ( New_P( i-1, Who_To_Dec ) + Decrement(i), zeros(1,N_Dec) );
        end
        
        % 5) Update the action's probability vector P
        %    FOR ACTIONS MORE LIKELY TO BE REWARDED
        
        if N_Inc ~= 0
            
            % Increment the selected probabilities and ensure they are bound to 1.
            % New_P( i-1, Who_To_Inc ) is the previous iteration:
            New_P( i, Who_To_Inc ) = min ( New_P( i-1, Who_To_Inc ) + Increment(i), ones(1,N_Inc) );
            
            % It can happen that the vector of probabilities sums to more than 1.
            % This is unphysical and should be avoided. The reason why this
            % happens is that the delta is too large, so that a certain
            % probability value is < 1 before update and >1 after update.
            % If this happens we need to take care of it (see below)
            
            Pre_Inc_Rem_Sum(i) =  sum(New_P( i, 1:Ai-1)) + sum( New_P( i, Ai+1:end) );
            

            % Check if the learning should stop
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
            
            
            % Protect against unphysicality of the vector of probabilities
            
            if Pre_Inc_Rem_Sum(i) > 1
                
                % Before updating the vector of probabilities the sum of all 
                % the elements except the one selected was above unity. We
                % calculate the excess and reduce the \Delta of that
                % amount.
                
                Excess(i) = abs(1- Pre_Inc_Rem_Sum(i) );
                
                New_Delta(i) = New_Delta(1) - Excess(i);
                
                % Ensure the increment amout is not zero
                Consistent_Increment(i) =  max( New_Delta(i)./K(i), 0 );
                
                % Increment the probabilities accordingly
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

        % 5) Update the action's probability vector P
        %    FOR CHOSEN ACTION
        
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
                
        % 6) Update the estimates of the probabilities of being rewarded
       
        New_D(i,:) =  W ./ Z; 

        % 7) Error Check
        who_is_negative = New_P(i,:) < 0;
        
        N_negatives = sum(who_is_negative);
        
        Sum_Probability = sum (New_P(i,:));

        if ( N_negatives ~= 0 ) 
            error(['The error occured at iteration number %d. We had %d negative probability(ies).',...
                ' The remaining sum was %d',i,N_negatives,Remaining_Sum(i)] );
        end    
        if (Sum_Probability > (1 + 10^-12) )
           error('The actions probabilities at iteration %d add up to %d which is greater than 1.', i,Sum_Probability);
        end

        i = i+1;
      
    end  %8) Repeat from 1 until a certai threshold is met
    
     % 9) Terminate the execution of the algorithm and assume learning complete.
    
end
