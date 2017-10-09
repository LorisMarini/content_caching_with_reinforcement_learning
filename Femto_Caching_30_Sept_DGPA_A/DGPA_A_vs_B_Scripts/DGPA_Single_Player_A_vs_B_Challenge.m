%{
-------------------- DGPA_Single_Player_Test -----------------------------


--------------------------   AUTHORSHIP  ----------------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

--------------------------   DESCRIPTION   --------------------------------

This script is designed as a testbed to compare the performance of the
learning algorithm DGPA with the DGPA_B (also referred to as CI Conditional Inaction). 

SST is the size of the statistical
sample and represents how many times the same learning is performed to then 
infer statistics on the learning speed and accuracy. S_Dim is the dimension 
of the search space and represents among how many actions the learner should 
learn from. The optimal action is changed everytime so that all actions in 
the search space are at least once optimal actions. 

Key functions: DGPA_Single_Player, DGPA_Single_Player_B


% ----------------------------   CODE   -----------------------------------
%}

clear all
close all

SST = 100;     % Statistical Sample Dimension
Sdim = 20;     % Dimension of the Search Space
Var = 0.8;     % Var is the environment variance (noise). Cannot be 0.

Real_Optimal_Actions = 1:1:Sdim;

% Number of real optimal actions tested
nroa = length( Real_Optimal_Actions);

% Pre allocate variables for speed
a = zeros(Sdim, nroa, SST );

Learning_Time          = a;
Number_of_Iterations   = a;
Correct_Learnings      = a;
Incorrect_Learnings    = a;
Learning_Time_B        = a;
Number_of_Iterations_B = a;
Correct_Learnings_B    = a;
Incorrect_Learnings_B  = a;

Errors = 0;
Errors_B = 0;

for s = 2:1:Sdim
    
    S = 1:1:s; 
    
    for oa = 1:1:s  
        
        % Version A - Standard DGPA
        for i = 1:1: SST
            try
                tic;
                % Key function 
                [ Learned_Action, Iteration_Number ] = DGPA_Single_Player( S ,Real_Optimal_Actions(oa),Var ,1 ,10);
                
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
        
        % Version B - Conditional Inaction DGPA
        for i = 1:1: SST
            try
                tic;
                % Key function 
                [ Learned_Action_B, Iteration_Number_B ] = DGPA_Single_Player_B( S ,Real_Optimal_Actions(oa),Var ,1 ,10);
                
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

% ----------------------- VERSION A -----------------------

Average_Correct_Learnings_PO = sum( Correct_Learnings,3)./SST;
Average_Wrong_Learnings_PO = sum(Incorrect_Learnings,3)./SST;
Average_Iterations_PO = sum(Number_of_Iterations,3)./SST;
Average_Learning_Time_PO = sum(Learning_Time,3)./SST;

Samples_Size = (sum(Average_Iterations_PO ~= 0,2));
Avg_N_Iter = sum(Average_Iterations_PO,2)./ Samples_Size;

Avg_Correct_Learnings = sum(Average_Correct_Learnings_PO,2)./ Samples_Size;
Avg_Wrong_Learnings = sum(Average_Wrong_Learnings_PO,2)./ Samples_Size;

Avg_Accuracy_Index = (Avg_Correct_Learnings ./ (Avg_Correct_Learnings + Avg_Wrong_Learnings));

% ----------------------- VERSION B -----------------------

Average_Correct_Learnings_PO_B = sum( Correct_Learnings_B,3)./SST;
Average_Wrong_Learnings_PO_B = sum(Incorrect_Learnings_B,3)./SST;
Average_Iterations_PO_B = sum(Number_of_Iterations_B,3)./SST;
Average_Learning_Time_PO_B = sum(Learning_Time_B,3)./SST;

Samples_Size_B = (sum(Average_Iterations_PO_B ~= 0,2));
Avg_N_Iter_B = sum(Average_Iterations_PO_B,2)./ Samples_Size_B;

Avg_Correct_Learnings_B = sum(Average_Correct_Learnings_PO_B,2)./ Samples_Size_B;
Avg_Wrong_Learnings_B = sum(Average_Wrong_Learnings_PO_B,2)./ Samples_Size_B;

Avg_Accuracy_Index_B = (Avg_Correct_Learnings_B ./ (Avg_Correct_Learnings_B + Avg_Wrong_Learnings_B));


% -------------------------- PLOT -----------------------

x_axis = [1:1:Sdim];

subplot(2,1,1), plot(x_axis ,Avg_Accuracy_Index);
title(['DGPA Average Accuracy Comparison - MLE - Normal Environment (Variance ' num2str(Var)  '), Statistical Sample Size = ' num2str(SST) '.']);
xlabel('Action Space Size');
ylabel('Average Accuracy');
xlim = [2,Sdim];
grid on;
hold on;
subplot(2,1,1), plot(x_axis ,Avg_Accuracy_Index_B, 'r');


subplot(2,1,2), plot(x_axis ,Avg_N_Iter);
title(['DGPA Convergence Comparison - MLE - Normal Environment (Variance ' num2str(Var)  '), Statistical Sample Size = ' num2str(SST) '.']);
xlabel('Action Space Size');
ylabel('Average Number of Iterations');
xlim = [2,Sdim];
grid on;
hold on;
subplot(2,1,2), plot(x_axis ,Avg_N_Iter_B,'r');

% Save the Workspace.
save(['DGPA_Convergence_Comparison_MLE_Normal Environment_Variance_' num2str(Var)  '_Statistical_Sample_Size_' num2str(SST)]);

