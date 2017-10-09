
function [ Learning, Iterations ] = INITIALIZE_Game_Of_DGPA( Network_Delays, M, Popularities, Reward_Type, Initialization_Number )
% Single Cell model of DGPA-BASED FEMTO CACHING

%% INPUT
% Resolution: = 1;
% Initialization_Number: = 10;
% Distances_Matrix: Matrix of distances.
% M: Caching capability of a single helper;
% Alpha: Free space attenuation factor;
%% OUTPUT
% Learning: THE INITIALISED LEARNERS (KxH)
%% Author: Loris Marini
% Version: 1.0.1
% Date: 09/09/2014

H = size(Network_Delays,2) - 1;            % Number of Helpers in the cell
N = size(Network_Delays,1);                % Number of users in the cell
F = H*M;                                   % Total number of files that can be cached.
S = 1:1:F;                                 % S: Space of Actions. F: How many files we can offload from the BS.
Sources_Degree = sum(Network_Delays < Inf, 1); % Number of users connected to each provider (Helper or BS)
Users_Degree = sum(Network_Delays < Inf, 2);   % Number of providers (Helper or BS) each user is connected to.

N_Ini = Initialization_Number;          % Initial #Iterations for estimation of D from each Learner
INI_Positive_Feedbacks = zeros(M,H);    % Initialise Environmental Feedback

Learning = Allocate_Learners( H,M,F );  % Learners' Variables Intialisation
Lesser_Selected = 0;                    % Is a variable to control the action that has been selected the least.
Zeros = Inf;                            % Is a variable to control the learners with empty D.
ITER = 1;                               % Iteration Number
Min_Weighted_Delay = Inf*ones(1,N);
Min_Average_Weighted_Delay = Inf*ones(1,H);
Average_Weighted_Delay = zeros(1,H);

while (Lesser_Selected < N_Ini || Zeros > 0)
    
   
    if ITER > 10000
       error('Initialisation Failed.');
    end
    
    %% INI Learners Select Files in Parallel (same time)
    
    [Available_Files, New_Learning] = Learners_Files_Selection( S, Learning );
    Learning = New_Learning;
    
    %% INI Feedbacks from the users

    INI_Rewards = zeros(M,H+1);   % Cumulative Rewards for all users
    INI_Penalties = zeros(M,H+1); % Cumulative Penalties for all users
    
    switch Reward_Type
        case 'Best_File_Based_Reward'
            for n = 1:1:N
                User_Selections = User_NCA_Selection( n, S, Available_Files, Network_Delays);
                Weighted_Delay = User_Weighted_Delay( User_Selections, Popularities );
                Delay_Performance(ITER,n) = Weighted_Delay;
                User_Delays = Network_Delays(n,:);
                [ Current_Rewards, N_Rewards, Current_Penalties, N_Penalties ] = Best_File_Based_Reward( User_Delays, User_Selections, Popularities );
                INI_Rewards = INI_Rewards + Current_Rewards;
                INI_Penalties = INI_Penalties + Current_Penalties;
            end
            
        case 'Weighted_Delay_Based_Reward'
            for n = 1:1:N
                User_Selections = User_NCA_Selection( n, S, Available_Files, Network_Delays);
                Weighted_Delay = User_Weighted_Delay( User_Selections, Popularities );
                Delay_Performance(ITER,n) = Weighted_Delay;
                User_Delays = Network_Delays(n,:);
                
                if (Weighted_Delay < Min_Weighted_Delay(n))
                    Min_Weighted_Delay(n) = Weighted_Delay;
                end
                if( ITER == 1)
                    Current_Minima = Inf;
                elseif (ITER > 1)
                    Current_Minima = Min_Weighted_Delay(n);
                end
                [ Current_Rewards, N_Rewards, Current_Penalties, N_Penalties ]...
                    = Weighted_Delay_Based_Reward( User_Delays, User_Selections, Weighted_Delay, Current_Minima , Popularities );
                
                INI_Rewards = INI_Rewards + Current_Rewards;
                INI_Penalties = INI_Penalties + Current_Penalties;
            end
        case 'Average_Weighted_Delay_Based_Reward'
            
            Tot_Rewards = zeros(M,H+1,N);
            Tot_Penalties = zeros(M,H+1,N);
            
            for n = 1:1:N
                User_Selections(n,:,:,:) = User_NCA_Selection( n, S, Available_Files, Network_Delays);
                This_User_Selections = squeeze(User_Selections(n,:,:,:));
                Weighted_Delay(n) = User_Weighted_Delay( This_User_Selections, Popularities );
                Delay_Performance(ITER,n) = Weighted_Delay(n);
                User_Delays = Network_Delays(n,:);
                 
                [ Current_Rewards, N_Rewards, Current_Penalties, N_Penalties ] = Best_File_Based_Reward( User_Delays,This_User_Selections, Popularities );
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
    
    %% INI Learners Determine the Environment Feedback Democratically
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
            if (max(Learning(k,j).D) == 0) % All the elements of D of (k,j) are zero. Until Zeros is ~= 0 we should keep initerating.
                Zeros = Zeros+1;
            end
            if (sum(Learning(k,j).D == 0) > MaxNZeros)
                MaxNZeros = sum(Learning(k,j).D == 0);
            end
        end
    end
    
    %disp(['INITIALISATION: Iteration ' num2str(ITER) ]);
    %disp(['The weighted delays are: ' num2str(Delay_Performance(ITER,:)) '.There are ' num2str(Zeros) ' learners with a zero D vector.']);
    %disp(['The maximum number of zeros in the D vectors is: ' num2str(MaxNZeros) '.'])
    ITER = ITER +1;
end
Iterations = ITER - 1;
%{
    disp('-------------------------------------------------------------------');
    disp('=================== INITIALISATION COMPLETE. ======================');
    disp('-------------------------------------------------------------------');
%}    
    %% OUTPUT Variable
    Learning = Learning;
    Network_Delays = Network_Delays;
    
end

