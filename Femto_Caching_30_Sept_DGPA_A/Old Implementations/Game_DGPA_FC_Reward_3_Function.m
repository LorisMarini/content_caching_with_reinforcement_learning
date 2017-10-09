function [  Minimum_Weighted_Delay ] = Game_DGPA_FC_Reward_3_Function(DistancesMatrix, InitializationNumber,Resolution )
% !!!!!!!!!!!!!!!!!!!!!!!!!!!

%  DO NOT USE IT IT'S WRONG.
%  DO NOT USE IT IT'S WRONG.
%  DO NOT USE IT IT'S WRONG.

% !!!!!!!!!!!!!!!!!!!!!!!!!!!

% Typical call:
    % Resolution = 1;
    % InitializationNumber = 10;
    % DistancesMatrix = Matrice delle distanze generata altrove.
    % Delays = Game_DGPA_FC_Reward_3(DistancesMatrix, 10,1)
    
% !!!!!!!!!!!!!!!!!!!!!!!!!!!

%  DO NOT USE IT IT'S WRONG.
%  DO NOT USE IT IT'S WRONG.
%  DO NOT USE IT IT'S WRONG.

% !!!!!!!!!!!!!!!!!!!!!!!!!!!

% Author: Loris Marini 
% Version: 1.0.1 
% Date: 09/09/2014

% DELAYS and FADING
% We create a matrix of delays that represents the status of the cell.The
% first column of 'Network_Delays' represents the delays of each user when
% downloading from helper 1. The second column from h2 and the last column
% represents the delays with the BS.
% The delays are computed assuming only free space loss. If a delay is
% 'Inf' means that there is not connection. The last column of 'Network_Delays' 
% is for the users delays with the BS and should not contain Inf. 
% We assume that each user has knowledge of its corresponding row of the
% 'Network_Delays' and therefore there is perfect CSI.

% FEEDBACK POLICY
% 1) For the file choosen by Learner k,j and available in a
% certain number of alternative places, each user should
% express a reward or a penalty. A reward shall be given to
% the selected source, while a penalty to the discarded ones except the BS.
% Each Action from Each Learner (k,j) needs to receive a Feedback from each of the N user.
% How many feedbacks should each helper receive? Each helper should
% receive a feedback from each user he is connected to, for every file
% he has, any time it appears across the network.

%{
%% Workspace Initialisation
clear all
close all
%}

 
%% Radio-Cell Configuration

%Distances = [12, 15, Inf, 26, 39; 27, 20, 26, Inf, 38; Inf, 20, 25, 30, 69; 19, Inf, 23, 16, 64];
%Distances = [12, Inf, 26, 39; 27, 26, Inf, 38; Inf, 25, 30, 69; 19, Inf, 16, 64];

Distances = DistancesMatrix;

H = size(Distances,2) - 1;            % Number of Helpers in the cell
N = size(Distances,1);                % Number of users in the cell
M = 3;                                % Caching Capability of a single helper
F = H*M;                              % Total number of files that can be cached.
N_Learners = F;                       % Total number of Learners
N_Ini = InitializationNumber;         % Initial #Iterations for estimation of D from each Learner 
Alpha = 3;                            % Free space attenuation factor.

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

%% Compute Delays in the Network:
% Network_Delays contains the delays experience by each user when connected to each helper.
% Network_Delays(n,j) = DELAY of user n when downloading from source j, (1< j < H+1 )
% Network_Delays(:,end) represents the DELAYS with the Base Station. 

Network_Delays = zeros(N, H + 1); 
Max_Capacity = log2( 1 + 1);
SNR = 1./(Distances.^Alpha); 
Capacities = log2( 1+ SNR);
Network_Delays = Capacities.^-1;

%% Determination of the Search Space 
S = [1:1:F];      % Space of Actions = How many files we can cache.
r = length(S);    % Number of actions to search from.

%% Determination of Popularities: ZIPF Distribution Law.

% We assume that S corresponds to the action's ranking. Action 1 is then
% the most popular file. Action F is the least popular file.
Zipf_Exp = 0.4;
Popularities = 1./(S.^Zipf_Exp)./(sum((1./S).^Zipf_Exp));


%% Network Topology

Sources_Degree = sum(Network_Delays < Inf, 1);
Users_Degree = sum(Network_Delays < Inf, 2);

%% Initialise Environmental Feedback

INI_Positive_Feedbacks = zeros(M,H);

%% Learners' Variables Intialisation

P = (1/r).*ones(1,r);  % Actions Probability Vector uniformally initialised
D = zeros(1,r);                 % Reward Probability Estimates
W = zeros(1,r);                 % Incidence of Rewards 
Z = zeros(1,r);                 % Actions Selection Incidence
Res = Resolution;               % Resolution Parameter
Delta = 1/(r.*Res);             % Resulution Step
Initial_P = (1/r).*ones(1,r);   % Initial value probability mass function
Initial_D = zeros(1,F);         % Initial estimates of the reward probability.

Learner = struct('P',Initial_P, 'Ai',0, 'Z',zeros(1,r), 'W',zeros(1,r), 'D',Initial_D);
 for k = 1:1:M
     for j = 1:1:H
         Learning(k,j) = Learner;
     end
 end
    
%% INITIALISATION

Lesser_Selected = 0;
Minimum_Weighted_Delay = Inf.*ones(1,N);
ITER = 1;
Min_Reward = 0;
Zeros = Inf;

while (Lesser_Selected < N_Ini || Zeros > 0)
    
    if ITER > 10000
        stop = 1;
    end
  
    Available_Files = zeros(M,H);
    
    %% Learners Select Files in Parallel (same time) 
    for j = 1:1:H
        for k = 1:1:M    
 
            if sum(Learning(k,j).P~=(1/r))>0
                error('During the initialisation actions should be selected uniformally.');
            end
            
            Action = randsrc(1,1,[ S; Learning(k,j).P ]);
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

   %% Cumulative Feedbacks from all the users
    Rewards = zeros(M,H+1);   % Cumulative Rewards for all users
    Penalties = zeros(M,H+1); % Cumulative Penalties for all users
    
    for n = 1:1:N
        Tot_Feedbacks_Received = zeros(M,H+1);
        N_Rewards = zeros(M,H+1);
        N_Penalties = zeros(M,H+1);
    
      %% USER's SELECTION (NCA)
        % User_Selections(k,j,f) is a collection of matrixes (M x H+1)
        % containing information about the NCA choices that user n
        % takes given the set of files (actions) taken by the
        % learners. For file 'f=1' for example, User_Selections(k,j,1):
        % User_Selections(k,j,1) = 0:     Discarded (will be penalised);
        % User_Selections(k,j,1) = Inf:   Non existing alternative;
        % User_Selections(k,j,1) = Delay: If user n has selected learner (k,j) to download file 'f'.
              
        User_Selections = zeros(M,H+1,F);
        
        for f = 1:1:F       
            
            Altern_From_Helpers = (Available_Files == f);    % Options for downloading file 'f': 
            Alternatives = Altern_From_Helpers;
            Alternatives(:,H+1) = true;                      % The BS has all files, always.
            Delays = Inf.*ones(1,H+1);
            Delays(sum(Alternatives,1)~=0) = ...             % Delays between the user n and the helpers that can provide file 'f'.
                Network_Delays( n, sum(Alternatives,1)~=0 );
                    
            [Min_Delay, S_ID_Selected] = min(Delays);        % S_ID_Selected = Helper or BS slected to download file f;
            L_ID_Selected = Alternatives(:,S_ID_Selected);   % L_ID_Selected = Learner within 'S_ID_Selected' that can provide file f;  
           
            % Calculate the selections for current user n on file f:
            
            for j = 1:1:H+1
                
                 Redundancy_Flag = 0;                             % This flag controls the penalties to redundanct actions. 
                 
                for k = 1:1:M
                    
                    if ( j == S_ID_Selected) 
                        
                        % The selected source can be a helper or the base
                        % station BS. We differentiate the two cases.
                        
                        if ( S_ID_Selected < H+1)                                   % We've selected the content from a helper.    
                            if(sum(L_ID_Selected) > 1)                              % In helper j there is a redundancy. 
                                    
                                if (Alternatives(k,j) == 1 && ~Redundancy_Flag)  
                                    User_Selections(k,j,f) = Delays(S_ID_Selected); % Report the delay with the selected learner.   
                                    Redundancy_Flag = 1;                            % All other alternatives in this helper will be penalised (redundant).   

                                elseif (Alternatives(k,j) == 1 && Redundancy_Flag) 
                                    User_Selections(k,j,f) = 0;                     % Penalise redunant selections.

                                elseif(Alternatives(k,j) == 0)
                                    User_Selections(k,j,f) = Inf;                   % All other learners are set to Infinity.
                                end     
                                
                            elseif(sum(L_ID_Selected) == 1)                         % In helper j there is NOT redundancy. 
                                
                                if (Alternatives(k,j) == 1)
                                    User_Selections(k,j,f) = Delays(S_ID_Selected);
                                elseif(Alternatives(k,j) == 0)
                                    User_Selections(k,j,f) = Inf; 
                                end       
                            end
                                
                        elseif( S_ID_Selected == H+1)                           % We've selected the content from THE BASE STATION.
                            
                            User_Selections(k,j,f) = Delays(S_ID_Selected);     % Report the delay with the selected learner.
                            
                        end        
                    
                  % If this source is not the selected one, we penalise (0) all
                  % the learners who can provide the alternatives and set
                  % to infinity all the others.
                  
                    elseif ( j ~= S_ID_Selected)                                % Penalise all the alternatives
                        
                        if (Alternatives(k,j) == 1)  
                            User_Selections(k,j,f) = 0;                         % Discard Alternative Reachable Learners.
                            
                        elseif(Alternatives(k,j) == 0)
                            User_Selections(k,j,f) = Inf;                       % Set to Infinity all others.
                        end
                    end      
                    
                end % For all files k
            end % For all files j
        end % For all files f     
        
        % OUTPUT: 
        % User_Selections
        
        
        % Safety Control: We don't want user n to reward more than one
        % learner for a given file. It has to choose.
        for f=1:1:F
            if ( sum(sum(     User_Selections(:,1:end-1,f) > 1  &  User_Selections(:,1:end-1,f) < Inf       )) > 1 )
                error(['User ' num2str(n) ' cannot reward multiple learners for file ' num2str(f) ' .']);
            end
        end
        
        
        %% USER's AVERAGE DELAY
        % User n can now determine its own AVERGAE delay when downloading 
        % the F files from the NC (Nearest Content) sources resulting from the Game.
        
        Weighted_Delay = 0;   % Average delay for user 'n'.
        N_Files_Selected = 0; % Error Check.
        
        % INPUT:
        % User_Selections
        
        for f = 1:1:F
            FLG = 1;
            
            for j = 1:1:H+1
                for k = 1:1:M
                    if ( j < H+1)
                        % When user 'n' has selected file f from a lerner in one helper
                        if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                            Weighted_Delay = Weighted_Delay + Popularities(f).* User_Selections(k,j,f);
                            N_Files_Selected = N_Files_Selected +1;
                        end
                    elseif ( j == H+1 && FLG)
                        % When user 'n' has selected file f from the Base Station
                        if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                            Weighted_Delay = Weighted_Delay + Popularities(f).* User_Selections(k,j,f);
                            N_Files_Selected = N_Files_Selected +1;
                            FLG = 0;
                        end
                    end
                end % For all Lerners (k)
            end % For all Sources (j)
        end % For all F files (f)   
      
        if (N_Files_Selected ~= r)
           error(['User ' num2str(n) ' did not select all the files in S when evaluating his weighted delay.']); 
        end
        
        % OUTPUT: 
        % Weighted_Delay   
        
      %% USER's FEEDBACK TO THE SOURCES  
        % If its AVERAGE delay is lower than the storical minimum:
        % It rewards the learnes selected and penalises those discarded.
        
        % If its AVERAGE delay is higher than the storical minimum:
        % It penalises both the lerners selected and those discarded.
       
        if (Weighted_Delay <= Minimum_Weighted_Delay(n) )
        
            Minimum_Weighted_Delay(n) = Weighted_Delay; 
            
            for f = 1:1:F  % For all F files (f)
                for j = 1:1:H+1 % For all Sources (j)
                    for k = 1:1:M % For all Lerners (k)
                        if (User_Selections(k,j,f) > 0 && User_Selections(k,j,f) ~= Inf)
                            % Reward the selected learners
                            Rewards(k,j) = Rewards(k,j) + Popularities(f);
                            N_Rewards(k,j) = N_Rewards(k,j) +1;
                        elseif (User_Selections(k,j,f) == 0)
                            % Penalise the discarded alternative learners.
                            Penalties(k,j) = Penalties(k,j) + Popularities(f);
                            N_Penalties(k,j) = N_Penalties(k,j) +1;
                        elseif (User_Selections(k,j,f) == Inf)
                            % Do Nothing.
                        end
                    end
                end
            end   
        else
            for f = 1:1:F  % For all F files (f)
                for j = 1:1:H+1 % For all Sources (j)
                    for k = 1:1:M % For all Lerners (k)
                        
                        if (User_Selections(k,j,f) == Inf)
                            % Do nothing.
                        else
                            Penalties(k,j) = Penalties(k,j) + Popularities(f);
                            N_Penalties(k,j) = N_Penalties(k,j) +1;
                        end
                    end
                end
            end
        end
  
        % Feedbacks error control across the learners (No BS)
        Feedbacks_Received = N_Rewards + N_Penalties;
        if ( sum(sum(Feedbacks_Received(:,1:end-1) > 1 )) > 0)
            error('There is a problem with the users decision');
        end
            
        Delay_Performance(ITER,:) = Minimum_Weighted_Delay;
        
    end % Users
     
     %% Learners Determine the Environment Feedback Democratically
      % ' Env_Feedback' contains the number of times a learner (k,j) get a
      % reward from the system. 
      for j = 1:1:H
        for k = 1:1:M          
            Curr_Action = Learning(k,j).Ai;  
            if (Rewards(k,j) >= Penalties(k,j))
                INI_Positive_Feedbacks(k,j) = INI_Positive_Feedbacks(k,j) + 1;
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action) + 1;
                Learning(k,j).D(Curr_Action)  = Learning(k,j).W(Curr_Action) / Learning(k,j).Z(Curr_Action) ;
                
            else
                INI_Positive_Feedbacks(k,j) = INI_Positive_Feedbacks(k,j);
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action);
                Learning(k,j).D(Curr_Action)  = Learning(k,j).D(Curr_Action) ;
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
      
      Zeros = 0;
      Min_Reward = Inf;
      MaxNZeros = 0;
      for j = 1:1:H
          for k = 1:1:M
              if (max(Learning(k,j).D) == 0)
                  %All the elements of D of (k,j) are zero.
                  % Until Zeros is ~= 0 we should keep initialising.
                  Zeros = Zeros+1;
              end
              
              if (sum(Learning(k,j).D == 0) > MaxNZeros)
                  MaxNZeros = sum(Learning(k,j).D == 0);
              end
          end
      end
      
      disp(['INITIALISATION: Iteration ' num2str(ITER) ' The min weighted delais are: ' num2str(Minimum_Weighted_Delay) '.There are ' num2str(Zeros) ' learners with a zero D vector.']);
      disp(['The maximum number of zeros in the D vectors is: ' num2str(MaxNZeros) '.'])
      ITER = ITER +1;
      
end % Iterations 
disp('-------------------------------------------------------------------');
disp('=================== INITIALISATION COMPLETE. ======================');
disp('-------------------------------------------------------------------');









%% Ready To Start The Game. Initialise Variables.

P_Threshold = 0.9;                            % Convergence Threshold.
Total_Time = 0;                               % Time Initialisation.
Conv_Actions = zeros(M,H);                    % Speed Initialisation.
Minimum_Weighted_Delay = Inf.*ones(1,N);      % Stability Initialisation.
GITER = 1; 
% Stability Initialisation.

%% Initialise Environmental Feedback

GAME_Positive_Feedbacks = zeros(M,H);

while ~Check_Game_Convergence( Learning, P_Threshold )
    
    GITER = GITER;                            % Game Iteration Number   
    tic;                                      % Iteration Timing
    Action = 0;                               % Stability Initialisation.
    Available_Files = zeros(M,H);             % Stability Initialisation.
    
    %% Parallel Files selection (Actions)
    
    for j = 1:1:H
        for k = 1:1:M
            Action = randsrc(1,1,[ S; Learning(k,j).P ]); 
            
            %{
            % Make sure that there are no redundancies in the cache of each helper.     
            % PROBLEM: 
            % It seems that even if we try not to have redundancies, some 
            % learners will converge to the same file. As a result the while 
            % loop here gets stuck. 
            
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
    
   %% Cumulative Feedbacks from all the users
   
    Rewards = zeros(M,H+1);   % Cumulative Rewards for all users
    Penalties = zeros(M,H+1); % Cumulative Penalties for all users
    
    for n = 1:1:N
        
        Tot_Feedbacks_Received = zeros(M,H+1);
        N_Rewards = zeros(M,H+1);
        N_Penalties = zeros(M,H+1);   
        
        
      %% USER's SELECTION (NCA)
        % User_Selections(k,j,f) is a collection of matrixes (M x H+1)
        % containing information about the NCA choices that user n
        % takes given the set of files (actions) taken by the
        % learners. For file 'f=1' for example, User_Selections(k,j,1):
        % User_Selections(k,j,1) = 0:     Discarded (will be penalised);
        % User_Selections(k,j,1) = Inf:   Non existing alternative;
        % User_Selections(k,j,1) = Delay: If user n has selected learner (k,j) to download file 'f'.
              
        User_Selections = zeros(M,H+1,F);
        
        for f = 1:1:F       
            
            Altern_From_Helpers = (Available_Files == f);    % Options for downloading file 'f': 
            Alternatives = Altern_From_Helpers;
            Alternatives(:,H+1) = true;                      % The BS has all files, always.
            Delays = Inf.*ones(1,H+1);
            Delays(sum(Alternatives,1)~=0) = ...             % Delays between the user n and the helpers that can provide file 'f'.
                Network_Delays( n, sum(Alternatives,1)~=0 );
                    
            [Min_Delay, S_ID_Selected] = min(Delays);        % S_ID_Selected = Helper or BS slected to download file f;
            L_ID_Selected = Alternatives(:,S_ID_Selected);   % L_ID_Selected = Learner within 'S_ID_Selected' that can provide file f;  
           
            % Calculate the selections for current user n on file f:
            
            for j = 1:1:H+1
                
                 Redundancy_Flag = 0;                             % This flag controls the penalties to redundanct actions. 
                 
                for k = 1:1:M
                    
                    if ( j == S_ID_Selected) 
                        
                        % The selected source can be a helper or the base
                        % station BS. We differentiate the two cases.
                        
                        if ( S_ID_Selected < H+1)                                   % We've selected the content from a helper.    
                            if(sum(L_ID_Selected) > 1)                              % In helper j there is a redundancy. 
                                    
                                if (Alternatives(k,j) == 1 && ~Redundancy_Flag)  
                                    User_Selections(k,j,f) = Delays(S_ID_Selected); % Report the delay with the selected learner.   
                                    Redundancy_Flag = 1;                            % All other alternatives in this helper will be penalised (redundant).   

                                elseif (Alternatives(k,j) == 1 && Redundancy_Flag) 
                                    User_Selections(k,j,f) = 0;                     % Penalise redunant selections.

                                elseif(Alternatives(k,j) == 0)
                                    User_Selections(k,j,f) = Inf;                   % All other learners are set to Infinity.
                                end     
                                
                            elseif(sum(L_ID_Selected) == 1)                         % In helper j there is NOT redundancy. 
                                
                                if (Alternatives(k,j) == 1)
                                    User_Selections(k,j,f) = Delays(S_ID_Selected);
                                elseif(Alternatives(k,j) == 0)
                                    User_Selections(k,j,f) = Inf; 
                                end       
                            end
                                
                        elseif( S_ID_Selected == H+1)                           % We've selected the content from THE BASE STATION.
                            
                            User_Selections(k,j,f) = Delays(S_ID_Selected);     % Report the delay with the selected learner.
                            
                        end        
                    
                  % If this source is not the selected one, we penalise (0) all
                  % the learners who can provide the alternatives and set
                  % to infinity all the others.
                  
                    elseif ( j ~= S_ID_Selected)                                % Penalise all the alternatives
                        
                        if (Alternatives(k,j) == 1)  
                            User_Selections(k,j,f) = 0;                         % Discard Alternative Reachable Learners.
                            
                        elseif(Alternatives(k,j) == 0)
                            User_Selections(k,j,f) = Inf;                       % Set to Infinity all others.
                        end
                    end      
                end
            end
        end % For all files f
        
        
      %{
        %% USER's SELECTION (NCA)
        % User_Selections(k,j,f) is a collection of matrixes (M x H+1)
        % containing information about the NCA choices that user n
        % takes given the set of files (actions) taken by the
        % learners. For file 'f=1' for example, User_Selections(k,j,1):
        % User_Selections(k,j,1) = 0:     Discarded (will be penalised);
        % User_Selections(k,j,1) = Inf:   Non existing alternative;
        % User_Selections(k,j,1) = Delay: If user n has selected learner (k,j) to download file 'f'.
              
        User_Selections = zeros(M,H+1,F);
        
        for f = 1:1:F       
            
            Altern_From_Helpers = (Available_Files == f);    % Options for downloading file 'f': 
            Alternatives = Altern_From_Helpers;
            Alternatives(:,H+1) = true;                      % The BS has all files, always.
            Delays = Inf.*ones(1,H+1);
            Delays(sum(Alternatives,1)~=0) = ...             % Delays between the user n and the helpers that can provide file 'f'.
                Network_Delays( n, sum(Alternatives,1)~=0 );
                    
            [Min_Delay, S_ID_Selected] = min(Delays);        % S_ID_Selected = Helper or BS slected to download file f;
            L_ID_Selected = Alternatives(:,S_ID_Selected);   % L_ID_Selected = Learner within 'S_ID_Selected' that can provide file f;  
            Redundancy_Flag = 0;                             % This flag controls the penalties to redundanct actions. 

            for j = 1:1:H+1
                for k = 1:1:M
                    
                    if ( j == S_ID_Selected)
                        
                        if ( S_ID_Selected < H+1 && sum(L_ID_Selected) > 1)  % We know that in helper j there is a redundancy.
                            if (Alternatives(k,j) == 1 && ~Redundancy_Flag)  
                                User_Selections(k,j,f) = Delays(S_ID_Selected); % Report the delay with the selected learner.   
                                Redundancy_Flag = 1;                            % Change the status of the flag.   
                                
                            elseif (Alternatives(k,j) == 1 && Redundancy_Flag) 
                                User_Selections(k,j,f) = 0;                     % Penalise redunant selections.
                                
                            elseif(Alternatives(k,j) == 0)
                                User_Selections(k,j,f) = Inf;                   % All other learners are set to Infinity.
                            end
                        else
                            if (Alternatives(k,j) == 1)                         %There is no redundancy in the selected helper j.
                                User_Selections(k,j,f) = Delays(S_ID_Selected); % Report the delay with the selected learner.  
                                
                            elseif(Alternatives(k,j) == 0)
                                User_Selections(k,j,f) = Inf;                   % All other learners are set to Infinity.
                            end
                        end                      
                    elseif ( j ~= S_ID_Selected)
                        
                        if (Alternatives(k,j) == 1)  
                            User_Selections(k,j,f) = 0;                         % Discard Alternative Reachable Learners.
                            
                        elseif(Alternatives(k,j) == 0)
                            User_Selections(k,j,f) = Inf;                       % Set to Infinity all others.
                        end
                    end      
                end
            end
        end % For all files f
        
        
        %}
        
        
      %% USER's AVERAGE DELAY
        % User n can now determine its own AVERGAE delay when downloading 
        % the F files from the NC (Nearest Content) sources resulting from the Game.      
        Weighted_Delay = 0;   % Average delay for user 'n'.
        N_Files_Selected = 0; % Error Check.
        
        for f = 1:1:F
            FLG = 1;
            
            for j = 1:1:H+1
                for k = 1:1:M
                    if ( j < H+1)
                        % When user 'n' has selected file f from a lerner in one helper
                        if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                            Weighted_Delay = Weighted_Delay + Popularities(f).* User_Selections(k,j,f);
                            N_Files_Selected = N_Files_Selected +1;
                        end
                    elseif ( j == H+1 && FLG)
                        % When user 'n' has selected file f from the Base Station
                        if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                            Weighted_Delay = Weighted_Delay + Popularities(f).* User_Selections(k,j,f);
                            N_Files_Selected = N_Files_Selected +1;
                            FLG = 0;
                        end
                    end
                end % For all Lerners (k)
            end % For all Sources (j)
        end % For all F files (f) 
        
        %{
         if (N_Files_Selected ~= F)
            error('ERROR G_234: Number of selections does not mach the number of files F.');
        end
        %}
          
  %% USER's FEEDBACK TO THE SOURCES  
        % If its AVERAGE delay is lower than the storical minimum:
        % It rewards the learnes selected and penalises those discarded.
        
        % If its AVERAGE delay is higher than the storical minimum:
        % It penalises both the lerners selected and those discarded.
       
        if (Weighted_Delay <= Minimum_Weighted_Delay(n) )
        
            Minimum_Weighted_Delay(n) = Weighted_Delay; 
            
            for f = 1:1:F  % For all F files (f)
                for j = 1:1:H+1 % For all Sources (j)
                    for k = 1:1:M % For all Lerners (k)
                        if (User_Selections(k,j,f) > 0 && User_Selections(k,j,f) ~= Inf)
                            % Reward the selected learners
                            Rewards(k,j) = Rewards(k,j) + Popularities(f);
                            N_Rewards(k,j) = N_Rewards(k,j) +1;
                        elseif (User_Selections(k,j,f) == 0)
                            % Penalise the discarded alternative learners.
                            Penalties(k,j) = Penalties(k,j) + Popularities(f);
                            N_Penalties(k,j) = N_Penalties(k,j) +1;
                        elseif (User_Selections(k,j,f) == Inf)
                            % Do Nothing.
                        end
                    end
                end
            end   
        else
            for f = 1:1:F  % For all F files (f)
                for j = 1:1:H+1 % For all Sources (j)
                    for k = 1:1:M % For all Lerners (k)
                        
                        if (User_Selections(k,j,f) == Inf)
                            % Do nothing.
                        else
                            Penalties(k,j) = Penalties(k,j) + Popularities(f);
                            N_Penalties(k,j) = N_Penalties(k,j) +1;
                        end
                    end
                end
            end
        end
  
        % Feedbacks error control across the learners (No BS)
        Feedbacks_Received = N_Rewards + N_Penalties;
        if ( sum(sum(Feedbacks_Received(:,1:end-1) > 1 )) > 0)
            error('There is a problem with the users decision');
        end
        Delay_Performance(GITER,:) = Minimum_Weighted_Delay;
        
    end % Users
    % Output1 = Rewards  (Matrix)
    % Output2 = Penalties  (Matrix)
    
    %% Learners Determine the Environment Feedback Democratically
      % 'Env_Feedback(k,j)'= 1  -->  Larner(k,j) Rewarded.
      % 'Env_Feedback(k,j)'= 0  -->  Learner(k,j) Penalised.
      for j = 1:1:H
        for k = 1:1:M          
            Curr_Action = Learning(k,j).Ai;  
            if (Rewards(k,j) >= Penalties(k,j))
                GAME_Positive_Feedbacks(k,j) = GAME_Positive_Feedbacks(k,j) + 1;
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action) + 1;
                Learning(k,j).D(Curr_Action)  = Learning(k,j).W(Curr_Action) / Learning(k,j).Z(Curr_Action) ;
                
            else
                GAME_Positive_Feedbacks(k,j) = GAME_Positive_Feedbacks(k,j);
                Learning(k,j).W(Curr_Action) = Learning(k,j).W(Curr_Action);
                Learning(k,j).D(Curr_Action)  = Learning(k,j).D(Curr_Action) ;
            end
        end
      end
      
  %% Probabilities Vectors Update. 
    % For each learner (k,j), how many actions where more likely to be
    % rewarded then the choosen action?

    Increment = zeros(M,H);
    Decrement = zeros(M,H);
    Excess = zeros(M,H);
    Penalty_Decrements = zeros(M,H);
    New_Delta = zeros(M,H);
    Consistent_Increment = zeros(M,H);
    
    for j = 1:1:H
        for k = 1:1:M
            
            if (Conv_Actions(k,j) == 1) % Do nothing, because learner (k,j) has already converged.
                Conv_Actions(k,j) = 1;
                
            elseif (Conv_Actions(k,j) == 0) % Learner (k,j) has not converged yet. Let's update its P & D.

                A = Learning(k,j).Ai;
                Curr_D = Learning(k,j).D;
                Direction_Vector = zeros(1,r);
                
                % 1) Determine which action should be incremented/decremented.            
                     % Direction_Vector = 1    INCREMENT
                     % Direction_Vector = 0    DECREMENT
                     % Direction_Vector = Inf  INACTION.

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


                % 2) Determine Increments and Decrements:
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
                
                % 4) Check Convergence
                if ( Check_Learner_Convergence( Learning(k,j).P, P_Threshold ) )     
                    % The Larner (k,j) has converged to its final decision,therefore, we can update its P vector:
                    [ LC, New_Prob_Vector, Action] = Check_Learner_Convergence( Learning(k,j).P, P_Threshold );
                    Conv_Actions(k,j) = Action;
                    Learning(k,j).P = New_Prob_Vector;   
                    disp(['Learner ' num2str(k) ' in helper ' num2str(j) ' has converged to file ' num2str(Action) '.']);
                    Conv_Actions
                end

                % 5) Decrease probabilities:  
                
                Learning(k,j).P( Who_To_Dec ) = max ( Learning(k,j).P( Who_To_Dec ) + Decrement(k,j), zeros(1,N_Dec(k,j))  );

                % 6) Adjust for Excess: (An excess can occur because we force probabilities to be >= 0)  
                     % If this happens we calculate the excess and we uniformally subtract it from the 
                     % probabilities that we have just increased.
                
                if (sum( Learning(k,j).P ) > 1)               
                    Excess(k,j) = abs(sum( Learning(k,j).P ) - 1 );
                    Penalty_Decrements(k,j) = Excess(k,j)/ K(k,j);
                    Learning(k,j).P( Who_To_Inc ) =  Learning(k,j).P( Who_To_Inc ) -   Penalty_Decrements(k,j);
                end

                % 7) Consistency Check:
                if (sum(  Learning(k,j).P ) > 1+10^-12 )
                    error('It didnt work. ERROR 101');
                end
                if (sum( Learning(k,j).P ) < 1-10^-12 )
                    error('It didnt work. ERROR 102');
                end
                if ( (sum(Learning(k,j).P < 0) ~= 0 ))
                    error(['The error occured at iteration number ' num2str(Iteration) '. We had '  num2str(N_negatives) ' negative probability(ies). The remaining sum was ' num2str(Sum_Probability)]);
                end

                % 8) UPDATE D.
                Learning(k,j).D = Learning(k,j).W ./ Learning(k,j).Z;   
            end       
        end % Lerner k
    end % Source j
    
    
    
    %% Iteration Timing
    Time = toc;
    Total_Time = Total_Time + Time;  
    Average_Time = Total_Time/GITER;
    %disp(['This is iteration number ' num2str(GITER) '. The Average Time per iteration is: ' num2str(Average_Time) '.'] );
   
    GITER = GITER + 1;% Game Iteration
    
end % While Loop. Iterations.

     
disp('-------------------------------------------------------------------');
disp('=================== Convergence COMPLETE. ======================');
disp('-------------------------------------------------------------------');


for n = 1:1:N
   disp(num2str( Minimum_Weighted_Delay(n)));
end

disp(['The Average Delay perceived by the users is: ' num2str(sum(Minimum_Weighted_Delay(n))./N) '.']);
disp(['The game with ' num2str(M*H) ' players (that is file memories) converged after ' num2str(GITER) ' iterations.']);















end

