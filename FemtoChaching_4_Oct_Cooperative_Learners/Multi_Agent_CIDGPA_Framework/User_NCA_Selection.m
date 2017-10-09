function [ User_Selections ] = User_NCA_Selection( User_Number, S, Available_Files, Network_Delays)

%% USER's SELECTION (NCA)
% User_Selections(k,j,f) is a collection of matrixes (M x H+1)
% containing information about the NCA choices that user n
% takes given the set of files (actions) taken by the
% learners. For file 'f=1' for example, User_Selections(k,j,1):
% User_Selections(k,j,1) = 0:     Discarded (will be penalised);
% User_Selections(k,j,1) = Inf:   Non existing alternative;
% User_Selections(k,j,1) = Delay: If user n has selected learner (k,j) to download file 'f'.

F = length(S);
NP = size(Available_Files,2) + 1; % Number of Providers (H+1)
NL = size(Available_Files,1);     % Number of Learners
n = User_Number;

User_Selections = zeros(NL,NP,F);

for f = 1:1:F
    
    Altern_From_Helpers = (Available_Files == f);    % Options for downloading file 'f':
    Alternatives = Altern_From_Helpers;
    Alternatives(:,NP) = true;                       % The BS has all files, always.
    Delays = Inf.*ones(1,NP);
    Delays(sum(Alternatives,1)~=0) = Network_Delays( n, sum(Alternatives,1)~=0 ); % Delays between the user n and the helpers that can provide file 'f'.
    [Min_Delay, S_ID_Selected] = min(Delays);                                     % S_ID_Selected = Helper or BS slected to download file f;
    L_ID_Selected = Alternatives(:,S_ID_Selected);                                % L_ID_Selected = Learner within 'S_ID_Selected' that can provide file f;
    
    % Calculate the selections for current user n on file f:
    
    for j = 1:1:NP       
        Redundancy_Flag = 0;                                                      % This flag controls the penalties to redundanct actions.       
        for k = 1:1:NL    
            if ( j == S_ID_Selected)
                % The selected source can be a helper or the base
                % station BS. We differentiate the two cases. 
                if ( S_ID_Selected < NP)                                    % We've selected the content from a helper.
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
                elseif( S_ID_Selected == NP)                                % We've selected the content from THE BASE STATION.           
                    User_Selections(k,j,f) = Delays(S_ID_Selected);         % Report the delay with the selected learner.     
                end
                
                % If this source is not the selected one, we penalise (0) all
                % the learners who can provide the alternatives and set
                % to infinity all the others.
                
            elseif ( j ~= S_ID_Selected)                 
                if ( Alternatives(k,j) == 0)                                % Set to Infinity those who didn't cache this file.
                    User_Selections(k,j,f) = Inf;                           % NOT AN OPTION                         
                elseif( Alternatives(k,j) == 1)
                    if( Delays(j) == Inf )                                  % Set to Infinity those you are not connected to.
                        User_Selections(k,j,f) = Inf;                       % NOT CONNECTED                       
                    end
                elseif( Alternatives(k,j) == 1)
                    if( Delays(j) ~= Inf )                                  % Set to Infinity those you are not connected to.
                        User_Selections(k,j,f) = 0;                         % Discard Alternative Reachable by the Learners.
                    end
                end
            end
            
        end % For all files k
    end % For all files j
end % For all files f

% Safety Control: We don't want user n to reward more than one
% learner for a given file. It has to choose.
for f=1:1:F
    if ( sum(sum(     User_Selections(:,1:end-1,f) > 1  &  User_Selections(:,1:end-1,f) < Inf       )) > 1 )
        error(['User ' num2str(n) ' cannot reward multiple learners for file ' num2str(f) ' .']);
    end
end
end

