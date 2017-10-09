function [ D,W,Z ] = MRC_Initialisation( Z,W,S,P,Opt_Act,Env_Variance,N_Ini)



% While the total number of actions selected is less than Init_Iterations

while  sum( Z < N_Ini ) ~= 0
    
    % Select an action at random
    Ai = randsrc(1,1,[S;P]);
    % Increment the frequency for that action
    Z(Ai) = Z(Ai) + 1;
    % Update the beleif in the optimal action (effec tof noise)
    Current_Optimum = Determine_Optimal_Action(Opt_Act, Env_Variance, S);
    
    % Ask a feedback to the environment and update the counter for the
    % rewards if the action gets a positive feedback.
    
    if (Environment(Ai,S,Current_Optimum))
        W(Ai) = W(Ai) + 1;
    else
        W(Ai) = W(Ai);
    end
end
% At the end of the initialization calculate the estimates on the
% reward probability as the ratio W/Z (which is proven to be maximally
% likely):
D = W./Z;

end

