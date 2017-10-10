function [ New_Learning, New_Positive_Feedbacks ] = Determine_Environment_Feedback( Learning, Rewards, Penalties, Positive_Feedbacks)

%% INI Learners Determine the Environment Feedback Democratically
% 'Env_Feedback(k,j)'= 1  -->  Larner(k,j) Rewarded.
% 'Env_Feedback(k,j)'= 0  -->  Learner(k,j) Penalised.

H = size(Learning,2);
M = size(Learning,1);


for j = 1:1:H
    for k = 1:1:M
        Curr_Action = Learning(k,j).Ai;
        Space = Learning(k,j).Search_Space;
        
        if ( Rewards(k,j) > Penalties(k,j) )
            Positive_Feedbacks(k,j) = Positive_Feedbacks(k,j) + 1;
            [rw, cl] = find( Space == Curr_Action );
            Learning(k,j).W( cl ) = Learning(k,j).W( cl ) + 1;
            Who_to_Divide = Learning(k,j).W ~= 0;
            Learning(k,j).D(Who_to_Divide) = Learning(k,j).W(Who_to_Divide)./ Learning(k,j).Z(Who_to_Divide);
        else
            % Do nothing.
        end
    end
end


New_Learning = Learning;
New_Positive_Feedbacks = Positive_Feedbacks;

end

