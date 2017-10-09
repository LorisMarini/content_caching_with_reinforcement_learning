function [ Learning, Iterations ] = Fast_Game_DGPA_INITIALISATION( Network_Delays, M, Popularities, Reward_Type, Initialization_Number  )

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
% Version: 1.0.2
% Date: 24/09/2014

H = size(Network_Delays,2) - 1;              % Number of Helpers in the cell
N = size(Network_Delays,1);                  % Number of users in the cell
F = H*M;                                     % Total number of files that can be cached.
S = 1:1:F;                                   % S: Space of Actions. F: How many files we can offload from the BS.
N_Ini = Initialization_Number;               % Initial #Iterations for estimation of D from each Learner
INI_Positive_Feedbacks = zeros(M,H);         % Initialise Environmental Feedback
Learning = Allocate_Learners( H,M,F );       % Learners' Variables Intialisation
Lesser_Selected = 0;                         % Is a variable to control the action that has been selected the least.
Zeros = Inf;                                 % Is a variable to control the learners with empty D.
ITER = 1;                                    % Iteration Number
Min_Weighted_Delay = Inf*ones(1,N);
Min_Average_Weighted_Delay = Inf*ones(1,H);
Average_Weighted_Delay = zeros(1,H);
Minima = Inf*ones(1,N);
            
while (Lesser_Selected < N_Ini || Zeros > 0)
    
    t_iteration_start = tic;
   
    if ITER > 10000
       error('Initialisation Failed.');
    end
    
    %% INI Learners Select Files in Parallel (same time)
    
    Available_Files = zeros(M,H);
    
    for j = 1:1: H
        for k = 1:1:M
            Action = randsrc(1,1,[ S; Learning(k,j).P ]);
            Available_Files(k,j) = Action;
            Learning(k,j).Ai = Action;
            Learning(k,j).Z(Action) = Learning(k,j).Z(Action) + 1;
        end
    end

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
         %{   
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
           %} 
            
        case 'Weighted_Delay_Based_Reward'
          t_feedbacks_start = tic;
            Rewards = zeros(M,H+1);
            Penalties = zeros(M,H+1);
            NCA = zeros(5, F+1);
            NCA(1,1:end-1) = 10^6*Popularities; % Set the Popularities

            for n=1:1:N                
                Reachable_Sources = repmat( Network_Delays(n,:) ~= Inf ,3,1);
                for f = 1:1:F       
                    Alternative_From_Helpers = (Available_Files == f);
                    Alternative_From_Providers = horzcat( Alternative_From_Helpers, ones(M,1)==1);
                    Sources_To_Consider = Alternative_From_Providers & Reachable_Sources;
                    [K_ids, J_ids] = find(Sources_To_Consider);
                    Competitors = horzcat(K_ids, J_ids);
                    
                    Delay_Quotes = Inf*ones(1,H+1);
                    Potential_Providers = sum(Sources_To_Consider,1) > 0;
                    Delay_Quotes(Potential_Providers) = Network_Delays(n, Potential_Providers );
                    [Min_Delay, J_Best_Provider] = min(Delay_Quotes);
                    NCA(2,f) = Min_Delay;                                % Set The minimum Delay
                    NCA(end,f) = NCA(1,f)/NCA(2,f);                  % Set The value of the reward
                    
                    Row_Best_Provider = find( Competitors(:,2) == J_Best_Provider, 1);
                    K_Best_Provider = Competitors( Row_Best_Provider, 1 );
                    NCA(end-1,f) = J_Best_Provider;
                    NCA(end-2,f) = K_Best_Provider;
                             
                    Learners_To_Penalise = Competitors;
                    Learners_To_Penalise(Row_Best_Provider,:) = [];
                    K_To_Penalise = Learners_To_Penalise(:,1);
                    J_To_Penalise = Learners_To_Penalise(:,2);
                    
                    Penalty_Weights = Delay_Quotes( sub2ind( [1,H+1], J_To_Penalise) );
                    Penalties( sub2ind([M,H+1],K_To_Penalise,J_To_Penalise) ) = Penalties( sub2ind([M,H+1],K_To_Penalise,J_To_Penalise) ) + (10^6*Popularities(f)./Penalty_Weights)';
                   
                end
                NCA(end,1:end-1) = NCA(1,1:end-1)/NCA(2,1:end-1);
                NCA(:,end) = sum( NCA(1,1:end-1).* NCA(2,1:end-1)); % Weighted Delay user n
                
                % Reward b
                if NCA(1,end) <=  Minima(n)
                    Minima(n) = NCA(1,end);
                    K_To_Reward = NCA(3,1:end-1);
                    J_To_Reward = NCA(4,1:end-1);
                    RW = NCA(end,1:end-1);
                    Rewards( sub2ind([M,H+1],K_To_Reward, J_To_Reward) ) = Rewards( sub2ind([M,H+1],K_To_Reward, J_To_Reward) ) + RW;
                end
           end
           t_feedbacks = toc(t_feedbacks_start);
           
        %{
        case 'Weighted_Delay_Based_Reward'
          t_feedbacks_start = tic;
            Rewards = zeros(M,H+1);
            Penalties = zeros(M,H+1);
            NCA = zeros(5, F+1, N);
            NCA(1,1:F,:) = repmat(10^6*Popularities,N,1)'; % Set the Popularities

            for n=1:1:N                
                Reachable_Sources = repmat( Network_Delays(n,:) ~= Inf ,3,1);
                for f = 1:1:F
                    
                    Alternative_From_Helpers = (Available_Files == f);
                    Alternative_From_Providers = horzcat( Alternative_From_Helpers, ones(M,1)==1);
                    Sources_To_Consider = Alternative_From_Providers & Reachable_Sources;
                    [K_ids, J_ids] = find(Sources_To_Consider);
                    Competitors = horzcat(K_ids, J_ids);
                    
                    Delay_Quotes = Inf*ones(1,H+1);
                    Potential_Providers = [sum(Sources_To_Consider,1) > 0];
                    Delay_Quotes(Potential_Providers) = Network_Delays(n, Potential_Providers );
                    [Min_Delay, J_Best_Provider] = min(Delay_Quotes);
                    NCA(2,f,n) = Min_Delay;                                % Set The minimum Delay
                    NCA(end,f,n) = NCA(1,f,n)/NCA(2,f,n);                  % Set The value of the reward
                    
                    Row_Best_Provider = find( Competitors(:,2) == J_Best_Provider, 1);
                    K_Best_Provider = Competitors( Row_Best_Provider, 1 );
                    NCA(end-1,f,n) = J_Best_Provider;
                    NCA(end-2,f,n) = K_Best_Provider;
                             
                    Learners_To_Penalise = Competitors;
                    Learners_To_Penalise(Row_Best_Provider,:) = [];
                    K_To_Penalise = Learners_To_Penalise(:,1);
                    J_To_Penalise = Learners_To_Penalise(:,2);
                    
                    Penalty_Weights = Delay_Quotes( sub2ind( [1,H+1], J_To_Penalise) );
                    Penalties( sub2ind([M,H+1],K_To_Penalise,J_To_Penalise) ) = Penalties( sub2ind([M,H+1],K_To_Penalise,J_To_Penalise) ) + (10^6*Popularities(f)./Penalty_Weights)';
                    
                    NCA(end,f,n) = NCA(1,f,n)/NCA(2,f,n);
                end
                NCA(:,end,n) = sum( NCA(1,1:end-1,n).* NCA(2,1:end-1,n)); % Weighted Delay user n
                
                % Reward b
                if NCA(1,end,n) <=  Minima(n)
                    Minima(n) = NCA(1,end,n);
                    K_To_Reward = NCA(3,1:end-1,n);
                    J_To_Reward = NCA(4,1:end-1,n);
                    RW = NCA(end,1:end-1,n);
                    Rewards( sub2ind([M,H+1],K_To_Reward, J_To_Reward) ) = Rewards( sub2ind([M,H+1],K_To_Reward, J_To_Reward) ) + RW;
                end
           end
           t_feedbacks = toc(t_feedbacks_start);
            %}
    end

    %% INI Learners Determine the Environment Feedback Democratically
    % 'Env_Feedback(k,j)'= 1  -->  Larner(k,j) Rewarded.
    % 'Env_Feedback(k,j)'= 0  -->  Learner(k,j) Penalised.
    
    for j = 1:1:H
        for k = 1:1:M
            Curr_Action = Learning(k,j).Ai;
            if (Rewards(k,j) > Penalties(k,j))
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
    
    disp(['INITIALISATION: Iteration ' num2str(ITER) ]);
    %disp(['The weighted delays are: ' num2str(Delay_Performance(ITER,:)) '.There are ' num2str(Zeros) ' learners with a zero D vector.']);
    %disp(['The maximum number of zeros in the D vectors is: ' num2str(MaxNZeros) '.'])
    ITER = ITER +1;
    
    t_iteration = toc(t_iteration_start);
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


