%{
-------------------- DGPA_Single_Player_Test -----------------------------


--------------------------   AUTHORSHIP  ----------------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

--------------------------   DESCRIPTION   --------------------------------

This script is designed as a testbed to monitor the performance of the
learning algorithm DGPA (single player). SST is the size of the statistical
sample and represents how many times the same learning is performed to then 
infer statistics on the learning speed and accuracy. S_Dim is the dimension 
of the search space and represents among how many actions the learner should 
learn from. The optimal action is changed everytime so that all actions in 
the search space are at least once optimal actions. 

Key function: DGPA_Single_Player


% ----------------------------   CODE   -----------------------------------
%}

clear all;
close all;

SST = 10;     % Statistical Sample Dimension
Sdim = 10;   % Dimension of the Search Space
Var = 1;      % Noise in the environment

Real_Optimal_Actions = 1:1:Sdim;
N_R_O = length( Real_Optimal_Actions);

Learning_Time = zeros(Sdim, N_R_O, SST );
Number_of_Iterations = zeros(Sdim, N_R_O, SST );
Correct_Learnings = zeros(Sdim, N_R_O, SST );
Incorrect_Learnings = zeros(Sdim, N_R_O,SST );

Errors = 0;

for s =1:1: Sdim
    S = 1:1:s; 
    
    % Set a different action as the optimal one, and see if the learning is
    % effective.
    
    for oa = 1:1: s  
        
        for i = 1:1: SST
            try
                tic;
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
          
    end
end

Average_Correct_Learnings_PO = sum( Correct_Learnings,3)./SST;
Average_Wrong_Learnings_PO = sum(Incorrect_Learnings,3)./SST;
Average_Iterations_PO = sum(Number_of_Iterations,3)./SST;
Average_Learning_Time_PO = sum(Learning_Time,3)./SST;

Samples_Size = (sum(Average_Iterations_PO ~= 0,2));
Avg_N_Iter = sum(Average_Iterations_PO,2)./ Samples_Size;

Avg_Correct_Learnings = sum(Average_Correct_Learnings_PO,2)./ Samples_Size;
Avg_Wrong_Learnings = sum(Average_Wrong_Learnings_PO,2)./ Samples_Size;

Avg_Accuracy_Index = (Avg_Correct_Learnings ./ (Avg_Correct_Learnings + Avg_Wrong_Learnings));

%%  Plot Data

close all;
figure;
plot([1:1:Sdim],Avg_Accuracy_Index);

fontsize = 14;
set(gcf, 'Color', 'w');
xlabel(  'Size of action Space',  'FontSize', fontsize );
ylabel( 'Accuracy Index', 'FontSize', fontsize );
title ({'Avg Accuracy Index vs Space Size'}, 'FontSize', fontsize );

figure;
plot([1:1:Sdim],Avg_N_Iter);

fontsize = 14;
set(gcf, 'Color', 'w');
xlabel(  'Size of action Space',  'FontSize', fontsize );
ylabel( 'I', 'FontSize', fontsize );
title ('Avg N Iterations vs Space Size', 'FontSize', fontsize );



