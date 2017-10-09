function [ Rewards, N_Rewards, Penalties, N_Penalties ] = Weighted_Delay_Based_Reward( Curr_User_Delays, User_Selections, Weighted_Delay, Min_Weighted_Delay, Popularities )

%% USER's FEEDBACK TO THE SOURCES
% If the users' AVERAGE delay is lower than the storical minimum:
% It rewards the learnes selected and penalises those discarded.
% If the users' AVERAGE delay is higher than the storical minimum:
% It penalises both the lerners selected and those discarded.
% It is assumed here that the selection is done according to the NCA
% principle (User_Selections).


NL = size(User_Selections,1); % Number of Learners
NP = size(User_Selections,2); % Number of Providers (H+1)
F = size(User_Selections,3);  % Number of Files
        
Rewards = zeros(NL,NP);
Penalties = zeros(NL,NP);
N_Rewards = zeros(NL,NP);
N_Penalties  = zeros(NL,NP);

for f = 1:1:F            % For all F files (f)
    for j = 1:1:NP       % For all Sources (j)
        for k = 1:1:NL   % For all Lerners (k)
            
             This_Helper_Delay = Curr_User_Delays(j);
             
            if (Weighted_Delay <= Min_Weighted_Delay )
                              
                if (User_Selections(k,j,f) > 0 && User_Selections(k,j,f) ~= Inf)
                    % Reward the selected learners
                    Rewards(k,j) = Rewards(k,j) + Popularities(f)/This_Helper_Delay;
                    N_Rewards(k,j) = N_Rewards(k,j) +1;
                elseif (User_Selections(k,j,f) == 0)
                    % Penalise the discarded alternative learners.
                    Penalties(k,j) = Penalties(k,j) + Popularities(f)/This_Helper_Delay;
                    N_Penalties(k,j) = N_Penalties(k,j) +1;
                elseif (User_Selections(k,j,f) == Inf)
                    % Do Nothing.
                end
            else
                if (User_Selections(k,j,f) == Inf)
                    % Do nothing.
                else
                    Penalties(k,j) = Penalties(k,j) + Popularities(f)/This_Helper_Delay;
                    N_Penalties(k,j) = N_Penalties(k,j) +1;
                end   
            end
            
        end 
    end
end

% Feedbacks error control across the learners (No BS)
Feedbacks_Received = N_Rewards + N_Penalties;
if ( sum(sum(Feedbacks_Received(:,1:end-1) > 1 )) > 0)
    error('There is a problem with the users decision. Too many feedbacks');
end

end
    