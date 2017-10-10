function [ Net_AWD_Best_Random ] = Performance_Best_Random_Placement( Results_Directory, ...
                                   Network_Delays, Gamma_ZipF, SST_Vector_Best_Rnd, M, H, N, S  )

%{
--------------------------   AUTHORSHIP  ---------------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

---------------------------   DESCRIPTION   -------------------------------

This script 

----------------------------- DEPENDENCIES --------------------------------

Files_Popularities()
Random_Placement()
Network_A_W_D()

-------------------------------- INPUT  ----------------------------------
           
Results_Directory
           The Directory where figures are saved.    
Network_Delays
           Matrix NxH+1 of delays between users and providers
Gamma_ZipF
           The range of parameters to try for the Zip-f probability ditrib.

SST_Vector_Best_Rnd
           

M
           Number of file memories in each helper.
H
           Number of helpers in the cell
N
           Number of users in the cell
S
           Library of files to cache

-------------------------------- OUTPUT  ---------------------------------- 

Net_AWD_Best_Random
          
% -----------------------------   CODE   ---------------------------------
%}

% Extract number of parameters to the Zip-f distribution                                  
N_Gammas = size(Gamma_ZipF,2);

Net_AWD_Best_Random = zeros(1,N_Gammas);

% SST_Vector_Best_Rnd = [23 44 78 ...]

if ( size(SST_Vector_Best_Rnd,1) ~= 1 || size(SST_Vector_Best_Rnd,2) ~= N_Gammas)
    error('This script is designed to work with a specific number of iterations per each gamma.');
end

T = 0;
Counter = 1;
Allocations_Left = sum(SST_Vector_Best_Rnd);

for g = 1:1:N_Gammas
    
    min_nawd = Inf;
    [ Popularities ] = Files_Popularities( S, Gamma_ZipF(g) ); 
    Iterations = SST_Vector_Best_Rnd(g);
    
    for i = 1:1:Iterations
        
        tstart = tic;
        [ ALLOC_Random ] = Random_Placement( M,H );
        this_net_awd = Network_A_W_D( Network_Delays, S, Popularities, ALLOC_Random );  
        
        if ( this_net_awd <= min_nawd)
            
            Net_AWD_Best_Random(g) = this_net_awd;
            min_nawd = this_net_awd;
        end
        
        T = T + toc(tstart);
        Avg_Time = T/Counter;
        Time_Left = Avg_Time*Allocations_Left;
        
        disp(['Best Random Placement. Time left: ', num2str(floor(Time_Left/60)), ...
             ' minutes and ' num2str(rem(Time_Left,60)) ' seconds.' ]);
         
        Counter = Counter + 1;
        Allocations_Left = Allocations_Left - 1;
    end
end

% Best out of the first 'SST_Vector_Best_Rnd' iterations:
Figure_Handle_1 = figure;
plot( Gamma_ZipF, Net_AWD_Best_Random, 'k'); 

Title = ['Best Placement out of iterations: ' num2str(SST_Vector_Best_Rnd)...
         ' H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
     
Fig_File_Name_1 = fullfile(Results_Directory, ['NAWD_Best_Random_N_' num2str(N) ...
                                               '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title);
xlabel('Gamma');
ylabel('NAWD');
grid on;
saveas( Figure_Handle_1,Fig_File_Name_1);

end



