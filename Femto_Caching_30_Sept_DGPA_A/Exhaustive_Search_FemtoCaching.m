function [ Optimal_Allocation ] = Exhaustive_Search_FemtoCaching( )
%{
-------------------------   AUTHORSHIP  -------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

-------------------------   DESCRIPTION   -------------------------

This function applies the exhaustive search to a particular cell with 4
users and 4 providers:

Distances = [12, Inf, 26, 39; 27, 26, Inf, 38; Inf, 25, 30, 69; 19, Inf, 16, 64];
Any attempt to make the network larger would result in insane processing
times. 

------------------------- INPUT PARAMETERS -------------------------

-- None --
 
------------------------- OUTPUT PARAMETERS -------------------------

-- Optimal_Allocation -- 
The optimal solution to this caching problem.

------------------------- EXAMPLE OF CALL -----------------------
%}


%% Network Topology:

Distances = [12, Inf, 26, 39; 27, 26, Inf, 38; Inf, 25, 30, 69; 19, Inf, 16, 64];

% Number of users in the network
N = size(Distances,1);

% Since the number of columns of Distances is the total number of content
% providers (Helpers + Base Station), the number of helpers is:
H = size(Distances,2) - 1;

% Cache size (that is number of files that each helper can cache)
M = 3;

% Total libary size F (That is the total number of files to cache)
F = H*M;

% Free space attenuation factor km^-1
Alpha = 3;          

% 
Available_Files = zeros(M,H);

% Sanity Check 1: ensure that Helpers are connected to more than one user so 
% that the caching problem is not trivial.

if min(sum(Distances < Inf,1)) <= 1
    error('Helpers must be connected to more than one user. Please use a valid Distances matrix.');
end

% Sanity Check 2: Helpers who can't help should be cut out.

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

Network_Delays = Network_Normalised_Delays(Distances, Alpha);


%% Determination of the Search Space 

S = [1:1:F];      % Space of Actions = one action epr each file in the library of size F.
r = length(S);    % Number of actions to search from.


%% Determination of Popularities: ZIPF Distribution Law.

% We assume that S is sorted in descending order by file popularity. 
% S(1) will be correspond to the action: "cache most poular file"; 
% S(2) corresponds to the action: "cache the second most popular file";
% S(3) ... etc. So we can calculate the corresponding popularities with a
% Zip-f distribution.

Zipf_Exp = 0.4;
Popularities = 1./(S.^Zipf_Exp)./(sum((1./S).^Zipf_Exp));

%% Exhaustive Search

% Pssible ways to take M files from the space of actions S:
N_Groups = nchoosek(S,M);

% Now each helper has N_Groups different choices. Since there are H
% helpers, the possible combinations of choices for all H helpers is:

N_Sets = unique(nchoosek(repmat(1:1:size(N_Groups,1), 1,H), H), 'rows');

% where repmat(1:1:size(N_Groups,1),1,H) represents all the actions
% available to the totality of H helpers and the unique function ensures
% that we only count file cachings that are different from each others.

% So if: N_Sets(400,:) is equal to 1,5,64 it means that the 400th netwrok
% configuration of H=3 helpers is such that the first helper take the first
% configuration in N_Groups, the seocnd takes the 5th and the third helper
% takes the 64th. 
     
     
if (size(N_Sets,1) ~= size(N_Groups,1)^H)
    error('mismatch');
end
if (size(N_Sets,2) ~= H)
    error('mismatch');
end


Min_Weighted_Delay = Inf*ones(1,N);

% Loop through all possible caching configurations and:

% 1. Extract the 'Available_Files' matrix containing the files that would be 
% cached by each helper in that setup 's'.

% 2. Simulate the user selections based on the criteria of nearest conent available. 
% If the same file is available in two locations the choie will be for the file 
% that can be downloaded faster (minimum latency).



for s = 1:1:size(N_Sets,1)
    
    % ----------------------------- 1 ------------------------------ 
    for h=1:1:H
        % Extract the cach configuration for helper h in set s:
        Current_Group = N_Sets(s,h); 
        % keep track of the files cached by each helper and therefore
        % available in the cell:
        Available_Files(:,h) = N_Groups(Current_Group,:);    
    end
    
    % ----------------------------  2 ------------------------------
    for n = 1:1:N
        
       % Apply user choices accoridng to the Nearest Content Available - NCA
       
       [ User_Selections ] = User_NCA_Selection( n, S, Available_Files, Network_Delays);
        
        %----  Calculate User average delay ----
        
        % Ok file complete and User_Selection populated. Now we can
        % calculate the average delay that the each user would experience
        % they wanted to sownload this file f:
         
        Weighted_Delay = zeros(1,N);   % Average delay for user 'n'.
        N_Files_Selected = 0; % Error Check.
        
        for f = 1:1:F
            flag = 1;
            
            for j = 1:1:H+1
                for k = 1:1:M
                    if ( j < H+1)
                        
                        % When user 'n' has selected file f from a lerner in one helper
                        if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                            
                            Weighted_Delay(n) = Weighted_Delay(n) + Popularities(f).* User_Selections(k,j,f);
                            
                            N_Files_Selected = N_Files_Selected +1;
                        end
                    elseif ( j == H+1 && flag)
                        
                        % When user 'n' has selected file f from the Base Station
                        if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                            
                            Weighted_Delay(n) = Weighted_Delay(n) + Popularities(f).* User_Selections(k,j,f);
                            N_Files_Selected = N_Files_Selected +1;
                            flag = 0;
                        end
                    end
                end % For all Lerners (k)
            end % For all Sources (j)
        end % For all F files (f) 
    
        if (Weighted_Delay(n) < Min_Weighted_Delay(n))
            
            % Update Weighted Delay if it is the smallest you have got so far
            Min_Weighted_Delay(n) = Weighted_Delay(n);
            
        end
    end  
    
    disp(['Set ' num2str(s) '/' num2str(size(N_Sets,1)) ' has been examined.',...
         ' The set of minimum weighted delays reached so far is: ' num2str(Min_Weighted_Delay) '.' ]);
end
    


end

