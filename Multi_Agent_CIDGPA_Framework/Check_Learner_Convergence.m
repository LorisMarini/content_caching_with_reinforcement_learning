function [ Learner_Convergence New_P_Vector Action] = Check_Learner_Convergence( Search_Space, Probability_Vector, P_Threshold )
% Check_Learner_Convergence
% This function is intended to check whether a Learner of the DGPA Game has
% converged or not after incrementing the probabilities. Note that the
% function handles the case of more than one action probability grater than
% one. In that case, the function will assume that the converged action is
% the one with the highest probability P.

% Learner_Convergence: Boolean Logic value. True = Learner has converged.
% False = The learner has not converged.

% New_P_Vector: Empty if Learner_Convergence = 0; However, if convergence 
% has been reached, New_P_Vector has probability 0 for all actions except
% the action to which the Lerner has converged (which has probability 1).


Learner_Convergence = 0;
Convergence_Index = Probability_Vector >= P_Threshold;
Converged_Actions = sum(Convergence_Index);

    if Converged_Actions == 0
           Learner_Convergence = 0;
           New_P_Vector = [];
           Action = NaN;
           
    elseif Converged_Actions == 1
           Learner_Convergence = true;
           New_P_Vector = Probability_Vector;
           New_P_Vector( Convergence_Index ) = 1;
           New_P_Vector( ~Convergence_Index ) = 0;
           [~, cl] = find(Convergence_Index);
           Action = Search_Space(cl);
           
           if sum(New_P_Vector) ~= 1
               error('The Function Check_Learner_Convergence is ill designed.');
           end           
    elseif Converged_Actions > 1
           Learner_Convergence = true;
           [rw, cl] = max(Probability_Vector(Convergence_Index));
           Action = Search_Space(cl);
           New_P_Vector = zeros(1,length(Probability_Vector));
           New_P_Vector(Action) = 1;
           if sum(New_P_Vector) ~= 1
               error('The Function Check_Learner_Convergence is ill designed.');
           end
    end
end

