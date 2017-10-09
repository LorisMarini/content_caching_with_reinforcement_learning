function [ Beta ] = Environment(Ai, S, Current_Opt )
% Environment is a function which determines a feedback Beta {0,1} based on 
% the Ai-th action selected from a set of action S, when the Current_Opt action is the optimal. 
% If Beta == 1 action Ai is rewarded, if Beta == 0 action Ai is penalised. 

%% Reward/Penalty Policy and Feedback Determination.

if Ai == Current_Opt   
    Beta = 1;   
else   
    Beta = 0; 
end


end

