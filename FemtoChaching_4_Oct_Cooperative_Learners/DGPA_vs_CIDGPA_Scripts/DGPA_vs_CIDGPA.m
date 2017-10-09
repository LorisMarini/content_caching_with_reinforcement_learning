% DGPA_Single_Player_Test
clear all
close all
% matlabpool('open',4);

SST = 100;     % Statistical Sample Dimension
S_Dim = 20;    % Dimension of the Search Space
Var = 0.8;   % Var is the environment variance (noise). Cannot be 0.

Real_Optimal_Actions = 1:1:S_Dim;
N_R_O = length( Real_Optimal_Actions);

Learning_Time = zeros(S_Dim, N_R_O, SST );
Number_of_Iterations = zeros(S_Dim, N_R_O, SST );
Correct_Learnings = zeros(S_Dim, N_R_O, SST );
Incorrect_Learnings = zeros(S_Dim, N_R_O,SST );

Learning_Time_B = zeros(S_Dim, N_R_O, SST );
Number_of_Iterations_B = zeros(S_Dim, N_R_O, SST );
Correct_Learnings_B = zeros(S_Dim, N_R_O, SST );
Incorrect_Learnings_B = zeros(S_Dim, N_R_O,SST );

Errors = 0;
Errors_B = 0;

for s = 2:1:S_Dim
    
    S = 1:1:s; 
    
    for oa = 1:1:s  
        
        % Version A
        for i = 1:1: SST
            try
                tic;
                [ Learned_Action, Iteration_Number ] = DGPA( S ,Real_Optimal_Actions(oa),Var ,1 ,10);
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
        
        % Version B
        for i = 1:1: SST
            try
                tic;
                [ Learned_Action_B, Iteration_Number_B ] = CI_DGPA( S ,Real_Optimal_Actions(oa),Var ,1 ,10);
                Learning_Time_B(s,oa,i) = toc;
                Number_of_Iterations_B(s,oa,i) = Iteration_Number_B;
                               
                if (Learned_Action_B == Real_Optimal_Actions(oa) )
                    Correct_Learnings_B(s,oa,i) = Correct_Learnings_B(s,oa,i) + 1;
                else
                    Incorrect_Learnings_B(s,oa,i) = Incorrect_Learnings_B(s,oa,i) + 1;
                end
                
            catch
                Errors_B = Errors_B + 1;
            end
        end 
        
    end
end

% matlabpool('close');

% VERSION A
Average_Correct_Learnings_PO = sum( Correct_Learnings,3)./SST;
Average_Wrong_Learnings_PO = sum(Incorrect_Learnings,3)./SST;
Average_Iterations_PO = sum(Number_of_Iterations,3)./SST;
Average_Learning_Time_PO = sum(Learning_Time,3)./SST;

Samples_Size = (sum(Average_Iterations_PO ~= 0,2));
Avg_N_Iter = sum(Average_Iterations_PO,2)./ Samples_Size;

Avg_Correct_Learnings = sum(Average_Correct_Learnings_PO,2)./ Samples_Size;
Avg_Wrong_Learnings = sum(Average_Wrong_Learnings_PO,2)./ Samples_Size;

Avg_Accuracy_Index = (Avg_Correct_Learnings ./ (Avg_Correct_Learnings + Avg_Wrong_Learnings));
x_axis = [1:1:S_Dim];

% VERSION B

Average_Correct_Learnings_PO_B = sum( Correct_Learnings_B,3)./SST;
Average_Wrong_Learnings_PO_B = sum(Incorrect_Learnings_B,3)./SST;
Average_Iterations_PO_B = sum(Number_of_Iterations_B,3)./SST;
Average_Learning_Time_PO_B = sum(Learning_Time_B,3)./SST;

Samples_Size_B = (sum(Average_Iterations_PO_B ~= 0,2));
Avg_N_Iter_B = sum(Average_Iterations_PO_B,2)./ Samples_Size_B;

Avg_Correct_Learnings_B = sum(Average_Correct_Learnings_PO_B,2)./ Samples_Size_B;
Avg_Wrong_Learnings_B = sum(Average_Wrong_Learnings_PO_B,2)./ Samples_Size_B;

Avg_Accuracy_Index_B = (Avg_Correct_Learnings_B ./ (Avg_Correct_Learnings_B + Avg_Wrong_Learnings_B));

subplot(2,1,1), plot(x_axis ,Avg_Accuracy_Index);
title(['DGPA Average Accuracy Comparison - MLE - Normal Environment (Variance ' num2str(Var)  '), Statistical Sample Size = ' num2str(SST) '.']);
xlabel('Action Space Size');
ylabel('Average Accuracy');
xlim = [2,S_Dim];
grid on;
hold on;
subplot(2,1,1), plot(x_axis ,Avg_Accuracy_Index_B, 'r');


subplot(2,1,2), plot(x_axis ,Avg_N_Iter);
title(['DGPA Convergence Comparison - MLE - Normal Environment (Variance ' num2str(Var)  '), Statistical Sample Size = ' num2str(SST) '.']);
xlabel('Action Space Size');
ylabel('Average Number of Iterations');
xlim = [2,S_Dim];
grid on;
hold on;
subplot(2,1,2), plot(x_axis ,Avg_N_Iter_B,'r');

% Save the Workspace.
save(['DGPA_Convergence_Comparison_MLE_Normal Environment_Variance_' num2str(Var)  '_Statistical_Sample_Size_' num2str(SST)]);

