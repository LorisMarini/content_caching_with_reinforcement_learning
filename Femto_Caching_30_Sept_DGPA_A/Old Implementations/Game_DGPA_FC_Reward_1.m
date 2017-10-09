% Single Cell model of DGPA-BASED FEMTO CACHING

% Author: Loris Marini 
% Version: 1.0.0 
% Date: 29/08/2014

% DELAYS and FADING
% We create a matrix of delays that represents the status of the cell.The
% first column of 'Network_Status' represents the delays of each user when
% downloading from helper 1. The second column from h2 and the last column
% represents the delays with the BS.
% The delays are computed assuming only free space loss. If a delay is
% 'Inf' means that there is not connection. The last column of 'Network_Status' 
% is for the users delays with the BS and should not contain Inf. 
% We assume that each user has knowledge of its corresponding row of the
% 'Network_Status' and therefore there is perfect CSI.

% FEEDBACK POLICY
% 1) For the file choosen by Learner k,j and available in a
% certain number of alternative places, each user should
% express a reward or a penalty. A reward shall be given to
% the selected source, while a penalty to the discarded ones except the BS.
% Each Action from Each Learner (k,j) needs to receive a Feedback from each of the N user.
% How many feedbacks should each helper receive? Each helper should
% receive a feedback from each user he is connected to, for every file
% he has, any time it appears across the network.

%% Workspace Initialisation
clear all
close all

%% Radio-Cell Configuration

H = 4;              % Number of Helpers in the cell
N = 4;              % Number of users in the cell
M = 5;              % Caching Capability of a single helper
F = H*M;            % Total number of files that can be cached.
N_Learners = M*H;   % Total number of Learners
N_Ini = 10;         % Initial #Iterations for estimation of D from each Learner 

Network_Status = zeros(N, H + 1); 

% FIRST SCENARIO:  Distances = [12, 15, Inf, 26, 39; 27, 20, 26, Inf, 18; Inf, 20, Inf, 19, 69; 19, Inf, Inf, 16, 64];
% SECOND SCENARIO: Distances = [12, 15, Inf, 26, 39; 27, 20, 17, Inf, 18; Inf, 20, Inf, 19, 69; 19, Inf, Inf, 16, 64];

% Distances = [12, 15, Inf, 26, 39; 27, 20, 26, Inf, 38; Inf, 20, 25, 19, 69; 19, Inf, 23, 16, 64];

Distances = [12, 15, Inf, 26, 39; 27, 20, 26, Inf, 38; Inf, 20, 25, 30, 69; 19, Inf, 23, 16, 64];

%% Check 1: Minimum Degree
if min(sum(Distances < Inf,1)) <= 1
    error('There are helpers connected only to a single user.');
end
%% Check 2: Helpers who can't help should be cut out.
for j=1:1:H
    for n=1:1:N
        if(Distances(n,j) > Distances(n,end) && Distances(n,j)< Inf)
            Distances(n,j) = Inf;
           disp(['In the network provided Helper ' num2str(j) ...
               'cannot be of any help to user ' num2str(j) ...
               '. The corresponding delay has been adjusted to Inf.']);
        end 
    end
end

Network_Status = Distances./90;
Env_Feedback = zeros(M,H);

%% Network Topology

Sources_Degree = sum(Network_Status < Inf, 1);
Users_Degree = sum(Network_Status < Inf, 2);

%% Determination of the Search Space 
S = [1:1:F];      % Space of Actions = How many files we can cache.
r = length(S);    % Number of actions to search from.

%% Determination of Popularities: ZIPF Distribution Law.

% We assume that S corresponds to the action's ranking. Action 1 is then
% the most popular file. Action F is the least popular file.
Zipf_Exp = 0.4;
Popularities = 1./(S.^Zipf_Exp)./(sum((1./S).^Zipf_Exp));

%% Learners Variables Intialisation

P = (1/r).*ones(1,r);  % Actions Probability Vector uniformally initialised
D = zeros(1,r);     % Reward Probability Estimates
W = zeros(1,r);     % Incidence of Rewards 
Z = zeros(1,r);     % Actions Selection Incidence
Res = 1;            % Resolution Parameter
Delta = 1/(r.*Res); % Resulution Step

Initial_P = (1/r).*ones(1,r);
Initial_D = zeros(1,F);

Learner = struct('P',Initial_P, 'Ai',0, 'Z',zeros(1,r), 'W',zeros(1,r), 'D',Initial_D);
 for k = 1:1:M
        for j = 1:1:H
            Learning(k,j) = Learner;
        end
 end
    
%% LEARNERS INITIALISATION

Lesser_Selected = 0;

while Lesser_Selected < N_Ini 
    
    %% Learners Select Files in Parallel (same time) 
    for j = 1:1:H
        for k = 1:1:M    
 Action = randsrc(1,1,[ S; Learning(k,j).P ]);
          
           % We allow redundancies in the cache of each helper.  
            Available_Files(k,j) = Action;
            Learning(k,j).Ai = Action;
            Learning(k,j).Z(Action) = Learning(k,j).Z(Action) + 1;
        end
    end
    
   %% Estimation of Number of Feedbacks Expected at each Helper
   
    N_Feedbacks = zeros(1,H);
    for j = 1:1:H
        for k = 1:1:M
            N_temp = sum(sum(Available_Files == Available_Files(k,j) ));
            N_Feedbacks(j) = N_Feedbacks(j) + N_temp;
        end
    end
    Tot_Feedbacks_Expected = Sources_Degree(1:end-1).*N_Feedbacks;
    
    % Reward and Penalty vectors must be initialised at every Iteration.
    Rewards = zeros(M, H+1);  
    Penalties = zeros(M, H+1);  
    
    % Variables to check consistency.
    N_RW = zeros(1,H+1);
    N_PEN = zeros(1,H+1);
    
    %% Learners Receive Feedbacks in Parallel (same time)

     for j = 1:1:H
        for k = 1:1:M    
            
            File_Popularity = Popularities(Learning(k,j).Ai); % Knowledge available to the Helpers.
            
            % Options for downloading file 'Learning(k,j).Ai' : ....
            
            Altern_From_Helpers = (Available_Files == Available_Files(k,j));
            Alternatives = Altern_From_Helpers;
            Alternatives(:,H+1) = true;
            
            % For file 'Learning(k,j).Ai', collect feedbacks form each user.
            % Rewards and penalties are proportional to the file popularity.
            % For each file, each user rewards the best learner (NCA)
            % and penalises the discarded learners that could have provided
            % that file.
            
            for n = 1:1:N
                
                % Connectivity check. User n must be connected to the
                % alternative sources:
                
                Delays = Inf.*ones(1,H+1);
                Delays(sum(Alternatives,1)~=0) = Network_Status( n, sum(Alternatives,1)~=0 ); 
                
                % Among the pool of sources that can provide the same file choosen by Learner (k,j), 
                % each user chooses according to the NCA (Nearest Content Available) policy.
                % If 'Sel_Source'= H+1 the users selected the BS.  
                
                [Min_Delay Sel_Source] = min(Delays);
                
                % LEARNER TO REWARD (ONLY ONE PER USER):   
                S_ID_To_Reward = Sel_Source;                 % Provider or S_ID (Either Helper of BS)
                L_ID_To_Reward = Alternatives(:,Sel_Source); % Learner within the S_ID to be rewarded.
                
                Rewards(L_ID_To_Reward, S_ID_To_Reward) = Rewards(L_ID_To_Reward, S_ID_To_Reward) + File_Popularity; 
                N_RW(S_ID_To_Reward) =  N_RW(S_ID_To_Reward) +sum(L_ID_To_Reward);
                
                %{
                if (Sel_Source < H+1 && sum(L_ID_To_Reward) > 1)  % Error check only for the helpers.
                    error('Attempt to reward more than one learner. Only one learner in each S_ID can be rewarded from a given user at a given time.');
                end  
                %}
                
                % LEARNERS TO PENALISE (CAN BE MORE THAN ONE PER USER):
                % Check vector Delays to see who is connected to user n, and penalise all learners 
                % who could offer the file but weren't selected by user n.
                for ss = 1:1:H+1      
                    if (Delays(ss) < Inf)
                        if (ss ~= Sel_Source)                                     
                            L_ID_To_Penalise = Alternatives(:,ss);
                            Penalties(L_ID_To_Penalise, ss) =  Penalties( L_ID_To_Penalise, ss) + File_Popularity; 
                            N_PEN(ss) =  N_PEN(ss) + sum(L_ID_To_Penalise); 
                        end
                    end
                end
                
            end % Number of Users
        end % Number of Learners
     end % Number of Helpers (1:H)
     
     % Coherence Check.
     Tot_Feedbacks_Received = N_RW + N_PEN;
     if ( sum(Tot_Feedbacks_Received(1:end-1) ~= Tot_Feedbacks_Expected) )
         error('There is a problem with the users decisionn');
     end
      
     %% Learners Determine the Environment Feedback Democratically
      % ' Env_Feedback' contains the number of times a learner (k,j) get a
      % reward from the system. 
      for j = 1:1:H
        for k = 1:1:M          
            Curr_Action = Learning(k,j).Ai;  
            if (Rewards(k,j) > Penalties(k,j))
                Env_Feedback(k,j) = Env_Feedback(k,j) + 1;
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action) + 1;
                Learning(k,j).D = Learning(k,j).W ./ Learning(k,j).Z;
            else
                Env_Feedback(k,j) = Env_Feedback(k,j);
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action);
                Learning(k,j).D = Learning(k,j).D;
            end
        end
      end
     
      Curr_min = min(Learning(1,1).Z);
      for j = 1:1:H
        for k = 1:1:M          
            if (min(Learning(k,j).Z) < Curr_min)
            Curr_min = min(Learning(k,j).Z);
            end
        end
      end
     
      Lesser_Selected = Curr_min;
      
end % TIME

% INITIALISATION COMPLETE. 

P_Threshold = 0.9;
Total_Time = 0;
Iteration = 1;

%% Ready To Start The Game. Initialise Variables.
Conv_Actions = zeros(M,H);

while ~Check_Game_Convergence( Learning, P_Threshold )
    
    tic;
    Action = 0;
    Available_Files = zeros(M,H);
    
    %% Parallel Files selection (Actions)
    for j = 1:1:H
        for k = 1:1:M
            Action = randsrc(1,1,[ S; Learning(k,j).P ]);
          
            % Make sure that there are no redundancies in the cache of each helper.  
            
            % PROBLEM: 
            % It seems that even if we try not to have redundancies, some 
            % learners will converge to the same file. As a result the while 
            % loop here gets stuck. 
            %{
            II = Conv_Actions(:,j)~=0;
            II(k) = false;       
            while  (  sum(Action == Conv_Actions(II,j)) > 0  ||  sum(Action == Available_Files(:,j)) > 0  )
                Action = randsrc(1,1,[ S; Learning(k,j).P ]);
                Available_Files
                Conv_Actions
                disp(['Learner ' num2str(k) ' in helper ' num2str(j) ' is trying to chose a suitable file...']);
                disp(['Current Selection: ' num2str(Action) '.']);                
            end
            %}
            Available_Files(k,j) = Action;
            Learning(k,j).Ai = Action;
            Learning(k,j).Z(Action) = Learning(k,j).Z(Action) + 1;
        end
    end
    
    %% Computation of Number of Feedbacks Expected
    N_Feedbacks = zeros(1,H);
    for j = 1:1:H
        for k = 1:1:M
            N_temp = sum(sum(Available_Files == Available_Files(k,j) ));
            N_Feedbacks(j) = N_Feedbacks(j) + N_temp;
        end
    end
    Tot_Feedbacks_Expected = Sources_Degree(1:end-1).*N_Feedbacks;
    % This vectors must be initialised at every Iteration.
    Rewards = zeros(M, H+1);
    Penalties = zeros(M, H+1);
    % Variables to check consistency.
    N_RW = zeros(1,H+1);
    N_PEN = zeros(1,H+1);
    
    %% Learners Receive Feedbacks in Parallel (same time)
    % We go through the action taken by each lerner.
    
    for j = 1:1:H
        for k = 1:1:M
            
            File_Popularity = Popularities(Learning(k,j).Ai);
            
            % Starting from the first file, users look if the same file is
            % available at any reachable helper.
            Altern_From_Helpers = (Available_Files == Available_Files(k,j));
            Alternatives = Altern_From_Helpers;
            Alternatives(:,H+1) = true;
            
            % For file 'Learning(k,j).Ai', collect feedbacks form each user.
            % Rewards and penalties are proportional to the file popularity.
            % For each file, each user rewards the best learner (NCA)
            % and penalises the discarded learners that could have provided
            % that file.
            
            for n = 1:1:N
                % Select only the users who can connect to the helper
                % who contains the alternative. Give a Reward to those whose delay is minimum.
                Delays = Inf.*ones(1,H+1);
                Delays(sum(Alternatives,1)~=0) = Network_Status( n, sum(Alternatives,1)~=0 );
                
                % Among the pool of sources that can provide the same file choosen by Learner (k,j),
                % each users chooses only the provider who can deliver the content with the lowest delay.
                % 'Sel_Source' tells which provider the users decided to download the file from. If 'Sel_Source'= H+1 the users selected the BS.
                
                [Min_Delay Sel_Source] = min(Delays);
                
                % LEARNER TO REWARD (ONLY ONE PER USER):
                S_ID_To_Reward = Sel_Source;                      % Provider or S_ID (Either Helper of BS)
                L_ID_To_Reward = Alternatives(:,Sel_Source);      % Learner withing the S_ID who needs to be rewarded.
                
                Rewards(L_ID_To_Reward, S_ID_To_Reward) = Rewards(L_ID_To_Reward, S_ID_To_Reward) + File_Popularity;
                N_RW(S_ID_To_Reward) =  N_RW(S_ID_To_Reward) + sum(L_ID_To_Reward);
                
                if (Sel_Source < H+1 && sum(L_ID_To_Reward) > 1)  % Error check only for the helpers.
%                    error('Attempt to reward more than one learner. Only one learner in each S_ID can be rewarded by a given user.');
                end
                
                % LEARNERS TO PENALISE (CAN BE MORE THAN ONE PER USER):
                % Check vector Delays to see who is connected to user n, and penalise all learners
                % who could offer the file but weren't selected by user n.
                for ss = 1:1:H+1
                    if (Delays(ss) < Inf)
                        if (ss ~= Sel_Source)
                            L_ID_To_Penalise = Alternatives(:,ss);
                            Penalties(L_ID_To_Penalise, ss) =  Penalties( L_ID_To_Penalise, ss) + File_Popularity;
                            N_PEN(ss) =  N_PEN(ss) + sum(L_ID_To_Penalise);
                            %N_PEN(ss) =  N_PEN(ss) + 1;
                        end
                    end
                end
                
            end % Number of Users
        end % Number of Learners
    end % Number of Helpers (1:H)
    
    % Coherence Check.
    Tot_Feedbacks_Received = N_RW + N_PEN;
    if ( sum(Tot_Feedbacks_Received(1:end-1) ~= Tot_Feedbacks_Expected) )
        error('There is a problem with the users decisionn');
    end
    
    %% Feedback Collection:
    % 'Env_Feedback' contains the number of times a learner (k,j) get a
    % reward from the system.Each learner chooses democratically.
    for j = 1:1:H
        for k = 1:1:M
            Curr_Action = Learning(k,j).Ai;
            if (Rewards(k,j) > Penalties(k,j))
                Env_Feedback(k,j) = Env_Feedback(k,j) + 1;
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action) + 1;
            else
                Env_Feedback(k,j) = Env_Feedback(k,j);
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action);
            end
        end
    end
    
    %% Probabilities Vectors Update.
    
    % For each learner (k,j), how many actions where more likely to be
    % rewarded then the choosen action?
    
    Pre_Inc_Rem_Sum = zeros(M,H);
    Post_Inc_Rem_Sum = zeros(M,H);
    Increment = zeros(M,H);
    Decrement = zeros(M,H);
    Excess = zeros(M,H);
    New_Delta = zeros(M,H);
    Consistent_Increment = zeros(M,H);
    
    for j = 1:1:H
        for k = 1:1:M
            
            % If learner (k,j) has not converged yet, update its P,D.
            
            if (Conv_Actions(k,j) == 0)
            
                A = Learning(k,j).Ai;
                Curr_D = Learning(k,j).D;

                % 1) Determine which action should be incremented and which decremented.
                Direction_Vector = zeros(1,r);
                % Vector element = 1  INCREMENT
                % Vector Element = 0  DECREMENT
                % Vector Element = Inf INACTION.

                Actions_Subset = Learning(k,j).P ~= 0;

                for o = 1:1:r
                    if (o ~= A)
                        if ( Curr_D(o) <= Curr_D(A)  &&  Learning(k,j).P( o )~= 0)
                            Direction_Vector(o) = 0;
                        elseif ( Curr_D(o) > Curr_D(A) &&  Learning(k,j).P( o )~= 0)                      
                            Direction_Vector(o) = 1;
                        elseif ( Learning(k,j).P( o ) == 0)  
                            Direction_Vector(o) = Inf;
                        end
                    elseif (o == A)
                        [Value, Position] = max( Curr_D(Actions_Subset) );
                        Number_of_Maxima = sum ( Curr_D(Actions_Subset) == Value );
                        if Number_of_Maxima == 1    % Single Maximum
                            if ( Curr_D(A) == Value )
                                Direction_Vector(o) = 1;
                            elseif ( Curr_D(A) ~= Value )
                                Direction_Vector(o) = 0;
                            end
                        elseif Number_of_Maxima > 1  % Multiple Maxima
                            Maxima =  Curr_D == Value;
                            Index = Maxima & Actions_Subset;
                            if sum ( Curr_D(Index) == Curr_D(A) ) > 1
                                Direction_Vector(o) = 1;
                            else
                                Direction_Vector(o) = 0;
                            end
                        end
                    end
                end
                K(k,j) = sum( Direction_Vector == 1 );
                N_Dec(k,j) = sum( Direction_Vector == 0 );
                Inactions(k,j) = sum(Direction_Vector == Inf);
                if (K(k,j) + N_Dec(k,j) +  Inactions(k,j)) ~= r
                    error('Inconsistency.');
                end 
                Who_To_Inc =  Direction_Vector == 1;
                Who_To_Dec =  Direction_Vector == 0;
                if (K(k,j) == 0)
                    error('There should always be an action to make more likely...');
                end


                % 2) Increments and Decrements.
                if ( K(k,j)~=0 )
                    Increment(k,j) = Delta./ K(k,j);
                else
                    Increment(k,j) = 0;
                end
                if ( N_Dec(k,j) ~=0)
                    Decrement(k,j) = - ( Delta ./ N_Dec(k,j) );
                else
                    Decrement(k,j) = 0;
                end

                % 3) Increase Probabilities          
                Learning(k,j).P( Who_To_Inc ) = Learning(k,j).P( Who_To_Inc ) + Increment(k,j);           
                if ( Check_Learner_Convergence( Learning(k,j).P, P_Threshold ) )     
                    % The Larner (k,j) has converged to its final decision.
                    % Therefore, we can update its P vector.
                    [ LC, New_Prob_Vector, Action] = Check_Learner_Convergence( Learning(k,j).P, P_Threshold );
                    Conv_Actions(k,j) = Action;
                    Learning(k,j).P = New_Prob_Vector;   
                    disp(['Learner ' num2str(k) ' in helper ' num2str(j) ' has converged to file ' num2str(Action) '.']);
                    Conv_Actions
                end

                % 4) Decrease probabilities:  
                Learning(k,j).P( Who_To_Dec ) = max ( Learning(k,j).P( Who_To_Dec ) + Decrement(k,j), zeros(1,N_Dec(k,j))  );

                % 5) Adjust for Excess:  
                if (sum( Learning(k,j).P ) > 1)               
                    Excess(k,j) = abs(sum( Learning(k,j).P ) - 1 );
                    Penalty_Decrements(k,j) = Excess(k,j)/ K(k,j);
                    Learning(k,j).P( Who_To_Inc ) =  Learning(k,j).P( Who_To_Inc ) -   Penalty_Decrements(k,j);
                end

               %% Consistency Check:
                if (sum(  Learning(k,j).P ) > 1+10^-12 )
                    error('It didnt work. ERROR 101');
                end
                if (sum( Learning(k,j).P ) < 1-10^-12 )
                    error('It didnt work. ERROR 102');
                end
                if ( (sum(Learning(k,j).P < 0) ~= 0 ))
                    error(['The error occured at iteration number ' num2str(Iteration) '. We had '  num2str(N_negatives) ' negative probability(ies). The remaining sum was ' num2str(Sum_Probability)]);
                end

               %% End of probability update. NOW UPDATE D
                Learning(k,j).D = Learning(k,j).W ./ Learning(k,j).Z;
                
            else
                % Do nothing, because learner (k,j) has already converged.
            end
                  
        end % Lerner k
    end % Source j
    
    
    
    %% Iteration Timing
    Time = toc;
    Total_Time = Total_Time + Time;  
    Average_Time = Total_Time/Iteration;
    disp(['This is iteration number ' num2str(Iteration) '. The Average Time per iteration is: ' num2str(Average_Time) '.'] );
    Iteration = Iteration + 1
    
end % While Loop. Iterations.

     
