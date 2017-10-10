function [ Available_Files, Helpers_Avg_Weighted_Delay ] = Greedy_Placement_Fast( Network_Delays, M, Popularities )

%{
-------------------------   AUTHORSHIP  -------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

-------------------------   DESCRIPTION   -------------------------

The algorithm implements a greedy caching strategy. It selects a random helper 
from the pool of available helpers and loops through all possible ways that
helper can cache M out of F files and picks the best one (the one that
minimizes the average latency that the users experience when requesting a 
file to that helper). It then goes on by selecting another helper at random
from the pool of remaning helpers until there are no more helpers left.
Random selections are always performed with unifrom probability mass
function. 

------------------------- INPUT PARAMETERS -------------------------

-- Network_Delays -- 
    The matrix of channel delays (H+1 columns, N rows)

-- M -- 
    Number of memory slots in each helper
       
-- Popularities --
    Array of file popularities.


------------------------- OUTPUT PARAMETERS -------------------------

-- Available_Files -- 
    The output of the caching problem, that is the files eventually cached
    by each helper.
    
-- helpers_avg_weighted_delay --
    the avg_weighted_delay for each helper after the optimization
    

------------------------- EXAMPLE OF CALL -----------------------

Alpha = 3;
Distances = [12, Inf, 26, 39; 27, 26, Inf, 38; Inf, 25, 30, 69; 19, Inf, 16, 64];
Network_Delays = Network_Normalised_Delays(Distances, Alpha);
M = 3;
Gamma_ZipF = 0.5;
F = 9; 
S = 1:1:F;
Popularities = Files_Popularities( S, Gamma_ZipF );
[ Available_Files, Helpers_Avg_Weighted_Delay ] = Greedy_Placement_Fast( Network_Delays, M, Popularities )

% ----------------------------   CODE     --------------------------
%}

% Number of Helpers.
H = size(Network_Delays,2) - 1;

% Number of Users.
N = size(Network_Delays,1);  

% Maximum number fo files that can be cached (also equal to the library
% size)
F = M*H;

% Space of actions (one action correspond to caching one file)
S = 1:1:F;         

% Initialize matrix of cached files
Available_Files = zeros(M,H);

% Number of possible ways we can select M actions (files) from the space of
% actions S:
N_Groups = nchoosek(S,M);

% 
Remaining_Helpers = 1:1:H;

% Initialize the mass probability function for the selection
P_Selection = (1/H)*ones(1,H);

% Initialize number of allocations left
Allocations_Left = size(N_Groups,1)*H;

% Prepare set of weights for the files (popularity)
Weights = repmat(Popularities, N, 1);

% Initialize helpers average weigthted delay
helpers_avg_weighted_delay = Inf*ones(1,H);


Counter = 1;
Total_Time = 0;

% Cycle through all helpers and apply greedy caching.

while size(Remaining_Helpers,2) ~= 0
    
    % Select first helper at random (uniform Probability Mass Function - PMF) 
    % from list of helpers:
    j = randsrc(1,1,[ Remaining_Helpers; P_Selection ]);
    
    % Eliminate the helper selected from the list of Helpers becasue it has just been selected
    Remaining_Helpers(Remaining_Helpers == j) = [];
    % Update the PMF to accout fro reduced number of helpers
    P_Selection = (1/size(Remaining_Helpers,2) )*ones(1, size(Remaining_Helpers,2) );
    
    % Initialize 'Min_Weighted_Delay'
    Min_Weighted_Delay = Inf;
    
    % For all possible ways in which the selected helper can cache M files,
    % take the users_NCA_delays for each user and for each file of the
    % library and calculate the global minimum of the avg_weighted_delay.
    % In doing so each helper thinks for itself and does not take into
    % accout what the other helpers do, from which the name "Greedy".
    
    for i = 1:1: size(N_Groups,1)
        
        tic; % Timing
        Available_Files(:,j) = N_Groups(i,:);
        users_NCA_delays = zeros(N,F);
        
        for f = 1:1:F
            % Collect download options for file f 
            options_from_helpers = (Available_Files == f); 
            
            who_hasnt_file_f = sum(options_from_helpers,1) == 0;
            
            Delays = Network_Delays;
            
            % Delays for helpers that cannot provide the file to infinite
            Delays(:, [who_hasnt_file_f, false] ) = Inf;
            
            % For each file in the library calculate the one with the
            % minimum latency for each user.
            users_NCA_delays(:,f) = min(Delays,[],2);
        end
        % Implement eq(4) of the original paper:
        % Marini, L., Li, J., & Li, Y. (n.d.). Distributed Caching based on Decentralized Learning Automata, 1–6.
        % There is one weighted delay for each user, and is calcuated as a
        % weighted sum of the users_NCA_delays.
        
        Weighted_Delays = sum(users_NCA_delays .* Weights,2);
        
        % Calculate how many users each helper is connected to:
        current_helper_degree = sum( Network_Delays(:,j) < Inf );
        
        % Calculate average weigthed delay for this helper.
        helpers_avg_weighted_delay(j) = sum( Weighted_Delays( Network_Delays(:,j) < Inf ) ) / current_helper_degree;
        
        % Update the minimum weighted delay if there is a need to
        if ( helpers_avg_weighted_delay(j) < Min_Weighted_Delay )
            
            Min_Weighted_Delay = helpers_avg_weighted_delay(j);
            
            % Also update what the best file allocation is.
            Current_Best_Allocation = N_Groups(i,:);
        end
        
        % Decrement Allocations left to exhamine
        Allocations_Left = Allocations_Left - 1;
        
        % Time monitoring
        Total_Time = Total_Time + toc;
        Avg_Partial_Time = Total_Time/Counter;
        Time_Left = Avg_Partial_Time * Allocations_Left;
        
        Counter = Counter + 1;
        disp(['GREEDY Placement. Time left: ' num2str(floor(Time_Left/60)),...
              ' minutes and ' num2str(rem(Time_Left,60)) ' seconds.']);
    end
    
    Available_Files(:,j) = Current_Best_Allocation;
    
end
end
