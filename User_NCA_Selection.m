function [ User_Selections ] = User_NCA_Selection( User_Number, S, Available_Files, Network_Delays)

%{
-------------------------   AUTHORSHIP  -------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

-------------------------   DESCRIPTION   -------------------------

Applies the Nearest Content Available - NCA policy at the user side in
order to understand where the user would download the files.

------------------------- INPUT PARAMETERS -------------------------

-- User Number --

    Integer representing the user identifier.

-- S -- 

    Action space. Is typically [1:1:F] where F is the size of the library (how
    many files can be cached).

-- Available_Files -- 

    A kxh matrix to keep track of the files cached by each helper and therefore
    available in the cell. The number of colums is equal to the number of
    helpers in the cell, and row equal to the number of memory slots that
    each helper has. 
        

-- Network_Delays --

    The matrix of latencies between each user and each provider.


------------------------- OUTPUT PARAMETERS -------------------------

-- User_Selections -- 

    User_Selections(k,j,f) is a collection of matrices, one for each
    file f. Since each provider caches a total of M files, and there
    are H+1 providers in the cell, these are (M x H+1) matrices. A
    column indicates the helper, and the row indicates which memory
    slot (max M). If user 'n' dhooses to download file 'f=1' from the
    k-th memory slot of provider 'j' then:

    User_Selections(k,j,1) = Delay:

    If instead it doesn't choose f=1 from that provider and slot,
    the User_selection for that file is set to zero so that the
    corresponding action may  be penalized:

    User_Selections(k,j,1) = 0:

    Finally if user 'n' cannot take file 'f' from this helper:

    User_Selections(k,j,1) = Inf:


------------------------- EXAMPLE OF CALL -----------------------


% ----------------------------   CODE     --------------------------
%}

                        
F = length(S);
NP = size(Available_Files,2) + 1; % Number of Providers (H+1)
NL = size(Available_Files,1);     % Number of Learners
n = User_Number;


User_Selections = zeros(NL,NP,F);


for f = 1:1:F
    
    % Retrive all the providers that can give this file (boolean)
    Altern_From_Helpers = (Available_Files == f);
    Alternatives = Altern_From_Helpers;
    
    % The BS has all files, always so set it to true
    Alternatives(:,NP) = true;
    
    % Initialize an array of delays for this user tryig to get this
    % file:
    Delays = Inf.*ones(1,NP);
    
    % Assignt the actual delay for all files available in (any) of the memory
    % slots of the providers:
    
    Delays(sum(Alternatives,1)~=0) = Network_Delays( n, sum(Alternatives,1)~=0 );
    
    % The file will is downloaded from the the provider with minimum delay:
    [Min_Delay, S_ID_Selected] = min(Delays);
    
    % L_ID_Selected = Learner within 'S_ID_Selected' that can provide file f;
    L_ID_Selected = Alternatives(:,S_ID_Selected);
            
 
    % -------------------------------------------------------------
    %                Populate User_Selections(k,j,f)
    % -------------------------------------------------------------
    
    % This flag controls the penalties to redundanct actions.
    Redundancy_Flag = 0;
    
    % For all providers and their memory slots:
    
    for j = 1:1:NP
        for k = 1:1:NL
            
            % If the provider is the one selected for this file f:
            if ( j == S_ID_Selected)
                
                % If there is a redundancy in helper j (same file cached twice)
                if ( S_ID_Selected < NP && sum(L_ID_Selected) > 1)
                    
                    % If the file is available and no redundancy has
                    % been fund yet:
                    if (Alternatives(k,j) == 1 && ~Redundancy_Flag)
                        
                        % Report the delay with the selected learner.
                        User_Selections(k,j,f) = Delays(S_ID_Selected);
                        % Change the status of the flag.
                        Redundancy_Flag = 1;
                        
                        % If file is available and is not the first
                        % time it is found redundant:
                        
                    elseif (Alternatives(k,j) == 1 && Redundancy_Flag)
                        % Penalise redundant selections.
                        User_Selections(k,j,f) = 0;
                        
                        % If file is NOT available:
                    elseif(Alternatives(k,j) == 0)
                        % All other learners are set to Infinity.
                        User_Selections(k,j,f) = Inf;
                    end
                    
                % Else if there is redundancy in the selected helper j:
                else
                    % If file is available
                    if (Alternatives(k,j) == 1)
                        % Report the delay with the selected learner.
                        User_Selections(k,j,f) = Delays(S_ID_Selected);
                        
                        % If file is NOT available
                    elseif(Alternatives(k,j) == 0)
                        % All other learners are set to Infinity.
                        User_Selections(k,j,f) = Inf;
                    end
                end
                
            elseif ( j ~= S_ID_Selected)
                
                if (Alternatives(k,j) == 1)
                    % Discard Alternative Reachable Learners.
                    User_Selections(k,j,f) = 0;
                    
                elseif(Alternatives(k,j) == 0)
                    % Set to Infinity all others.
                    User_Selections(k,j,f) = Inf;
                end
            end
            
        end % For all memory slots
    end % For all providers
end % For all files f



% Sanity Check: We don't want user n to reward more than one learner for a given file. It has to choose.

for f=1:1:F
    if ( sum(sum(     User_Selections(:,1:end-1,f) > 1  &  User_Selections(:,1:end-1,f) < Inf       )) > 1 )
        error(['User ' num2str(n) ' cannot reward multiple learners for file ' num2str(f) ' .']);
    end
end
end

