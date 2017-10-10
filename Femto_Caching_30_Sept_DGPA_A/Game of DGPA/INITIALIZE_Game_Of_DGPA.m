
function [ Learning, Iterations ] = INITIALIZE_Game_Of_DGPA( Network_Delays, M, Popularities, Reward_Type, Initialization_Number )
                                                         
%{
--------------------------   AUTHORSHIP  ---------------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

---------------------------   DESCRIPTION   -------------------------------

This script initializes a game of learning automata using DGPA
reinforcement learning.

----------------------------- DEPENDENCIES --------------------------------

Allocate_Learners(...)
Learners_Files_Selection(...)
User_NCA_Selection(...)
User_Weighted_Delay(...)
                                           
Best_File_Based_Reward()
Weighted_Delay_Based_Reward()

-------------------------------- INPUT  ----------------------------------
                                                
Network_Delays 
    Matrix of delays NxH+1 where N= number of users and H+1= Number of providers                                                         
M 
    Caching capability of a single helper;
                                                             
Popularities
    File Popularities
                                                             
Reward_Type
                                                                                                  
Initialization_Number
    See NITIALIZE_Game_Of_DGPA

-------------------------------- OUTPUT  ----------------------------------

Learning
    KxH Initialized learners

Iterations
    The actual number of iteratiosn that it took to initialize.

% -----------------------------   CODE   ---------------------------------
%}


H = size(Network_Delays,2) - 1;             % Number of Helpers in the cell
N = size(Network_Delays,1);                 % Number of users in the cell
F = H*M;                                    % Total number of files that can be cached.
S = 1:1:F;                                  % S: Space of Actions. F: How many files we can offload from the BS.

N_Ini = Initialization_Number;              % Initial #Iterations for estimation of D from each Learner
INI_Positive_Feedbacks = zeros(M,H);        % Initialise Environmental Feedback
Zeros = Inf;                                % Is a variable to control the learners with empty D.
ITER = 1;                                   % Iteration Number

Min_Weighted_Delay = Inf*ones(1,N);
Min_Average_Weighted_Delay = Inf*ones(1,H);
Average_Weighted_Delay = zeros(1,H);


% Loop until you make sure that all actions are selected at least N_Ini
% times. We do it with 'Lesser_Selected' which is a variable that controls 
% the action that has been selected the least.

Lesser_Selected = 0;   

% Initialize learners
Learning = Allocate_Learners( H,M,F );

while (Lesser_Selected < N_Ini || Zeros > 0)
    
    % If we reach 10^4 iterations we declared the initialization to have
    % failed. since actions are selected at random, for small numbers of
    % N_Ini and smalle sets convergence should be reached much earlier.
    % However, larger sets and N_Init might need to adjust this threshold
    % to a larger value.
    
    if ITER > 10000
       error('Initialisation Failed.');
    end
    
    % INI Learners Select Files in Parallel (same time)
    
    [Available_Files, New_Learning] = Learners_Files_Selection( S, Learning );
    Learning = New_Learning;
    
    % INI Feedbacks from the users:
    
    % Pre-allocate cumulative Rewards for all users
    % Pre-allocate cumulative Penalsties for all users
    INI_Rewards = zeros(M,H+1);   
    INI_Penalties = zeros(M,H+1); 
    
    switch Reward_Type
        
        case 'Best_File_Based_Reward'
                     
            for n = 1:1:N
                % Let each user choose which content they prefer using the
                % policy if nearest available (Nearest Content Available
                % NCA)
                user_selections = User_NCA_Selection( n, S, Available_Files, Network_Delays);
                
                % Calculate the weighted latency that the user would experience
                delays(ITER,n) = User_Weighted_Delay( user_selections, Popularities );
                
                this_user_delays = Network_Delays(n,:);
                
                % Let the user provide its own feedback
                [Current_Rewards, ~, Current_Penalties, ~ ] = Best_File_Based_Reward( this_user_delays, user_selections, Popularities );
                
                % Update Rewards and Penalties
                INI_Rewards = INI_Rewards + Current_Rewards;
                INI_Penalties = INI_Penalties + Current_Penalties;
                
            end
            
        case 'Weighted_Delay_Based_Reward'
            
            for n = 1:1:N
                
                % Let each user choose which content they prefer using the
                % policy if nearest available (Nearest Content Available
                % NCA)
                user_selections = User_NCA_Selection( n, S, Available_Files, Network_Delays);
                
                % Calculate the weighted latency that the user would experience
                delays(ITER,n) = User_Weighted_Delay( user_selections, Popularities );
                
                % Extract the actual delays that interest this user
                this_user_delays = Network_Delays(n,:);
                
                if (Weighted_Delay < Min_Weighted_Delay(n))
                    
                    Min_Weighted_Delay(n) = Weighted_Delay;
                end
                if( ITER == 1)
                    Current_Minima = Inf;
                elseif (ITER > 1)
                    Current_Minima = Min_Weighted_Delay(n);
                end
                
                [ Current_Rewards, ~, Current_Penalties, ~ ]...
                    = Weighted_Delay_Based_Reward( this_user_delays, user_selections, Weighted_Delay, Current_Minima , Popularities );
                
                % Update Rewards and Penalties
                INI_Rewards = INI_Rewards + Current_Rewards;
                INI_Penalties = INI_Penalties + Current_Penalties;
            end
        case 'Average_Weighted_Delay_Based_Reward'
            
            Tot_Rewards = zeros(M,H+1,N);
            Tot_Penalties = zeros(M,H+1,N);
            
            for n = 1:1:N
                % Let each user choose which content they prefer using the
                % policy if nearest available (Nearest Content Available
                % NCA)
                user_selections(n,:,:,:) = User_NCA_Selection( n, S, Available_Files, Network_Delays);
                This_User_Selections = squeeze(user_selections(n,:,:,:));
                
                delays(ITER,n)  = User_Weighted_Delay( This_User_Selections, Popularities );
                this_user_delays = Network_Delays(n,:);
                 
                [ Current_Rewards, ~, Current_Penalties, ~ ] = Best_File_Based_Reward( this_user_delays,This_User_Selections, Popularities );
                Tot_Rewards(:,:,n) = Current_Rewards;
                Tot_Penalties(:,:,n) = Current_Penalties;
            end
            
            for j=1:1:H
                % Take only the users connected to the learner (j,k)
                Who_to_Average = Network_Delays(:,j) ~= Inf;
                Average_Weighted_Delay(j) = sum( Weighted_Delay(Who_to_Average)) / sum(Who_to_Average);
                
                if (Average_Weighted_Delay(j) <= Min_Average_Weighted_Delay(j) )
                    
                    Min_Average_Weighted_Delay(j) = Average_Weighted_Delay(j);
                    INI_Rewards(:,j) = INI_Rewards(:,j) + sum( squeeze( Tot_Rewards(:,j,:) ),2); 
                    INI_Penalties(:,j) = INI_Penalties(:,j) + sum( squeeze( Tot_Penalties(:,j,:) ),2); 
                else
                    INI_Rewards(:,j) = INI_Rewards(:,j); 
                    INI_Penalties(:,j) = INI_Penalties(:,j) + sum( squeeze( Tot_Rewards(:,j,:) ),2) + sum( squeeze( Tot_Penalties(:,j,:) ),2); 
                end
            end
    end
    
    % INI Learners Determine the Environment Feedback Democratically:
    
    % 'Env_Feedback(k,j)'= 1  -->  Larner(k,j) Rewarded.
    % 'Env_Feedback(k,j)'= 0  -->  Learner(k,j) Penalised.
    for j = 1:1:H
        for k = 1:1:M
            Curr_Action = Learning(k,j).Ai;
            
            if (INI_Rewards(k,j) > INI_Penalties(k,j))
                
                INI_Positive_Feedbacks(k,j) = INI_Positive_Feedbacks(k,j) + 1;
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action) + 1;
                Who_to_Divide = Learning(k,j).W ~= 0;
                Learning(k,j).D(Who_to_Divide) = Learning(k,j).W(Who_to_Divide)./ Learning(k,j).Z(Who_to_Divide);
              
            else
                % Do nothing.
            end
        end
    end
    
    Curr_min = Inf;
    for j = 1:1:H
        for k = 1:1:M
            if (min(Learning(k,j).Z) < Curr_min)
                Curr_min = min(Learning(k,j).Z);
            end
        end
    end
    Lesser_Selected = Curr_min;
    
    Zeros = 0;
    MaxNZeros = 0;
    for j = 1:1:H
        for k = 1:1:M
            % All the elements of D of (k,j) are zero. Until Zeros is ~= 0 we should keep initerating.
            if (max(Learning(k,j).D) == 0)
                Zeros = Zeros+1;
            end
            if (sum(Learning(k,j).D == 0) > MaxNZeros)
                MaxNZeros = sum(Learning(k,j).D == 0);
            end
        end
    end
    
    disp(['INITIALISATION: Iteration ' num2str(ITER) ]);
    disp(['The weighted delays are: ' num2str(delays(ITER,:)) '.There are ' num2str(Zeros) ' learners with a zero D vector.']);
    disp(['The maximum number of zeros in the D vectors is: ' num2str(MaxNZeros) '.'])
    ITER = ITER +1;
end
Iterations = ITER - 1;

% OUTPUT Variables

Learning = Learning;
Network_Delays = Network_Delays;
    
end

