function [ OUT_Learning, OUT_Converged_Actions ] = New_P_Vectors_DGPA_B( Conv_Actions, IN_Learning, Delta, P_Threshold )


%% Probabilities Vectors Update.
% For each learner (k,j), how many actions where more likely to be
% rewarded then the choosen action?

M = size(IN_Learning,1);
H = size(IN_Learning,2);
r = size(IN_Learning(1,1).P,2);
K = zeros(M,H);
N_Dec = zeros(M,H);
Inactions = zeros(M,H);
Increment = zeros(M,H);
Decrement = zeros(M,H);
Excess = zeros(M,H);
Penalty_Decrements = zeros(M,H);

for j = 1:1:H
    for k = 1:1:M
        
        if (Conv_Actions(k,j) == 1) % Do nothing, because learner (k,j) has already converged.
            Conv_Actions(k,j) = 1;
            
        elseif (Conv_Actions(k,j) == 0) % Learner (k,j) has not converged yet. Let's update its P & D.
            
            A = IN_Learning(k,j).Ai;
            Curr_D = IN_Learning(k,j).D;
            Direction_Vector = zeros(1,r);
            
            % 1) Determine which action should be incremented/decremented.
            % Direction_Vector = 1    INCREMENT
            % Direction_Vector = 0    DECREMENT
            % Direction_Vector = Inf  INACTION.
            
            Actions_Subset = IN_Learning(k,j).P ~= 0;
            
            for o = 1:1:r
                if (o ~= A)
                    if ( Curr_D(o) <= Curr_D(A)  &&  IN_Learning(k,j).P( o )~= 0)
                        Direction_Vector(o) = 0;
                    elseif ( Curr_D(o) > Curr_D(A) &&  IN_Learning(k,j).P( o )~= 0)
                        Direction_Vector(o) = 1;
                    elseif ( IN_Learning(k,j).P( o ) == 0)
                        Direction_Vector(o) = Inf;
                    end
                elseif (o == A)
                    [Value, Position] = max( Curr_D(Actions_Subset) );
                    Number_of_Maxima = sum ( Curr_D(Actions_Subset) == Value );
                    if Number_of_Maxima == 1    % Single Maximum
                        if ( Curr_D(A) == Value )
                            Direction_Vector(o) = 1;
                        elseif ( Curr_D(A) ~= Value )
                            Direction_Vector(o) = 0;
                        end
                    elseif Number_of_Maxima > 1  % Multiple Maxima
                        Maxima =  Curr_D == Value;
                        Index = Maxima & Actions_Subset;
                        if sum ( Curr_D(Index) == Curr_D(A) ) > 1
                            Direction_Vector(o) = 1;
                        else
                            Direction_Vector(o) = 0;
                        end
                    end
                end
            end
            K(k,j) = sum( Direction_Vector == 1 );
            N_Dec(k,j) = sum( Direction_Vector == 0 );
            Inactions(k,j) = sum(Direction_Vector == Inf);
            if (K(k,j) + N_Dec(k,j) +  Inactions(k,j)) ~= r
                error('Inconsistency.');
            end
            Who_To_Inc =  Direction_Vector == 1;
            Who_To_Dec =  Direction_Vector == 0;
            if (K(k,j) == 0)
                error('There should always be an action to make more likely...');
            end
            
            
            % 2) Determine Increments and Decrements:
            if ( K(k,j)~=0 )
                Increment(k,j) = Delta./ K(k,j);
            else
                Increment(k,j) = 0;
            end
            if ( N_Dec(k,j) ~=0)
                Decrement(k,j) = - ( Delta ./ N_Dec(k,j) );
            else
                Decrement(k,j) = 0;
            end
            
            % 3) Increase Probabilities
            
            IN_Learning(k,j).P( Who_To_Inc ) = IN_Learning(k,j).P( Who_To_Inc ) + Increment(k,j);
            
            % 4) Check Convergence
            if ( Check_Learner_Convergence( IN_Learning(k,j).P, P_Threshold ) )
                % The Larner (k,j) has converged to its final decision,therefore, we can update its P vector:
                [ LC, New_Prob_Vector, Action] = Check_Learner_Convergence( IN_Learning(k,j).P, P_Threshold );
                Conv_Actions(k,j) = Action;
                IN_Learning(k,j).P = New_Prob_Vector;
                %disp(['Learner ' num2str(k) ' in helper ' num2str(j) ' has converged to file ' num2str(Action) '.']);
                %disp( num2str(Conv_Actions));
            end
            
            % 5) Decrease probabilities:
            
            IN_Learning(k,j).P( Who_To_Dec ) = max ( IN_Learning(k,j).P( Who_To_Dec ) + Decrement(k,j), zeros(1,N_Dec(k,j))  );
            
            % 6) Adjust for Excess: (An excess can occur because we force probabilities to be >= 0)
            % If this happens we calculate the excess and we uniformally subtract it from the
            % probabilities that we have just increased.
            
            if (sum( IN_Learning(k,j).P ) > 1)
                Excess(k,j) = abs(sum( IN_Learning(k,j).P ) - 1 );
                Penalty_Decrements(k,j) = Excess(k,j)/ K(k,j);
                IN_Learning(k,j).P( Who_To_Inc ) =  IN_Learning(k,j).P( Who_To_Inc ) -   Penalty_Decrements(k,j);
            end
            
            % 7) Consistency Check:
            if (sum(  IN_Learning(k,j).P ) > 1+10^-12 )
                error('It didnt work. ERROR 101');
            end
            if (sum( IN_Learning(k,j).P ) < 1-10^-12 )
                error('It didnt work. ERROR 102');
            end
            if ( (sum(IN_Learning(k,j).P < 0) ~= 0 ))
                error(['The error occured at iteration number ' num2str(Iteration) '. We had '  num2str(N_negatives) ' negative probability(ies). The remaining sum was ' num2str(Sum_Probability)]);
            end
            
            % 8) UPDATE D.
            IN_Learning(k,j).D = IN_Learning(k,j).W ./ IN_Learning(k,j).Z;
        end
    end % Lerner k
end % Source j
    
%% Assign output variables.

OUT_Learning = IN_Learning;
OUT_Converged_Actions = Conv_Actions;  

end

