function [ Optimal_Allocation ] = Exhaustive_Search_FemtoCaching( )
% Exhaustive search File Association.

%   This code runs only for this network set up:
%   Distances = [12, Inf, 26, 39; 27, 26, Inf, 38; Inf, 25, 30, 69; 19, Inf, 16, 64];

%% Network Topology:

Distances = [12, Inf, 26, 39; 27, 26, Inf, 38; Inf, 25, 30, 69; 19, Inf, 16, 64];
N = size(Distances,1);
H = size(Distances,2) - 1;
M = 3;
F = H*M;
Alpha = 3;          % Free space attenuation factor.
Available_Files = zeros(M,H);

% Check 1: Minimum Degree
if min(sum(Distances < Inf,1)) <= 1
    error('There are helpers connected only to a single user.');
end
% Check 2: Helpers who can't help should be cut out.
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

%% EXHAUSTIVE SEARCH

N_Groups = nchoosek(1:1:F,M);
N_Sets = unique(nchoosek(repmat(1:1:size(N_Groups,1), 1,H), H), 'rows');

if (size(N_Sets,1) ~= size(N_Groups,1)^H)
    error('mismatch');
end
if (size(N_Sets,2) ~= H)
    error('mismatch');
end


Min_Weighted_Delay = Inf*ones(1,N);


for s = 1:1:size(N_Sets,1)
    
    for h=1:1:H
        Current_Group = N_Sets(s,h); 
        Available_Files(:,h) = N_Groups(Current_Group,:);  % Column 1 is helper 1.  
    end
    
    for n = 1:1:N
    
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
        
      %% USER's AVERAGE DELAY
        % User n can now determine its own AVERGAE delay when downloading 
        % the F files from the NC (Nearest Content) sources resulting from the Game.
        
        Weighted_Delay = zeros(1,N);   % Average delay for user 'n'.
        N_Files_Selected = 0; % Error Check.
        
        for f = 1:1:F
            FLG = 1;
            
            for j = 1:1:H+1
                for k = 1:1:M
                    if ( j < H+1)
                        % When user 'n' has selected file f from a lerner in one helper
                        if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                            Weighted_Delay(n) = Weighted_Delay(n) + Popularities(f).* User_Selections(k,j,f);
                            N_Files_Selected = N_Files_Selected +1;
                        end
                    elseif ( j == H+1 && FLG)
                        % When user 'n' has selected file f from the Base Station
                        if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                            Weighted_Delay(n) = Weighted_Delay(n) + Popularities(f).* User_Selections(k,j,f);
                            N_Files_Selected = N_Files_Selected +1;
                            FLG = 0;
                        end
                    end
                end % For all Lerners (k)
            end % For all Sources (j)
        end % For all F files (f) 
    
        if (Weighted_Delay(n) < Min_Weighted_Delay(n))
            
            Min_Weighted_Delay(n) = Weighted_Delay(n);
            
        end
    end  
    
    disp(['Set ' num2str(s) '/' num2str(size(N_Sets,1)) ' has been examined. The set of minimum weighted delays reached so far is: ' num2str(Min_Weighted_Delay) '.' ]);
end
    


end

