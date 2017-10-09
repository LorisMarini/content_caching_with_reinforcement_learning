% DGPA_Single_Player_Test
clear all
close all
% matlabpool('open',4);

SST = 10;     % Statistical Sample Dimension
S_Dim = 10;   % Dimension of the Search Space

Real_Optimal_Actions = 1:1:S_Dim;
N_R_O = length( Real_Optimal_Actions);

Learning_Time = zeros(S_Dim, N_R_O, SST );
Number_of_Iterations = zeros(S_Dim, N_R_O, SST );
Correct_Learnings = zeros(S_Dim, N_R_O, SST );
Incorrect_Learnings = zeros(S_Dim, N_R_O,SST );

Errors = 0;

for s =1:1: S_Dim
    S = 1:1:s; 
    
    for oa = 1:1: s  
        
        for i = 1:1: SST
            try
                tic;
                [ Learned_Action, Iteration_Number ] = DGPA_Single_Player( S ,Real_Optimal_Actions(oa),1 ,1 ,10);
                Learning_Time(s,oa,i) = toc;
                Number_of_Iterations(s,oa,i) = Iteration_Number;

                if (Learned_Action == Real_Optimal_Actions(oa) )
                    Correct_Learnings(s,oa,i) = Correct_Learnings(s,oa,i) + 1;
                else
                    Incorrect_Learnings(s,oa,i) = Incorrect_Learnings(s,oa,i) + 1;
                end
                
            catch
                Errors = Errors + 1;       
            end 
        end
          
    end
end

% matlabpool('close');

Average_Correct_Learnings_PO = sum( Correct_Learnings,3)./SST;
Average_Wrong_Learnings_PO = sum(Incorrect_Learnings,3)./SST;
Average_Iterations_PO = sum(Number_of_Iterations,3)./SST;
Average_Learning_Time_PO = sum(Learning_Time,3)./SST;

Samples_Size = (sum(Average_Iterations_PO ~= 0,2));
Avg_N_Iter = sum(Average_Iterations_PO,2)./ Samples_Size;

Avg_Correct_Learnings = sum(Average_Correct_Learnings_PO,2)./ Samples_Size;
Avg_Wrong_Learnings = sum(Average_Wrong_Learnings_PO,2)./ Samples_Size;

Avg_Accuracy_Index = (Avg_Correct_Learnings ./ (Avg_Correct_Learnings + Avg_Wrong_Learnings));

plot([1:1:S_Dim],Avg_Accuracy_Index);
figure
plot([1:1:S_Dim],Avg_N_Iter);
