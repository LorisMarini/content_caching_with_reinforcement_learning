function [ Beta ] = Environment(a, S, Current_Opt )

% Environment is a function which determines the feedback Beta {0,1} based on 
% the a-th action selected from a set of action S. Current_Opt action is the optimal. 
% If Beta == 1 action 'a' is rewarded, if Beta == 0 action 'a' is penalised. 

%% Reward/Penalty Policy and Feedback Determination.

if a == Current_Opt   
    Beta = 1;   
else   
    Beta = 0; 
end


end

