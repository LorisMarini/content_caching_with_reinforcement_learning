function [ Failures, Avg_Net_AWD_DGPA, Avg_Iter_GAME_DGPA, Avg_Iter_INIT_DGPA, Iterations_DGPA_INIT, Iterations_DGPA_GAME ] ...
    = Performance_Learning( Learning_Setup, Results_Directory, Reward_Type, Network_Delays, Gamma_ZipF, SST_L, M, H, N, S, Ini_Number, Resolution, Conv_Prob_Th )

%{
--------------------------   AUTHORSHIP  ---------------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

---------------------------   DESCRIPTION   -------------------------------

This Function sets out to calculate the Average of the Average Weighted
Delay of the users for a set of Gammas (Zip-f distribution parameter) and
a statistical sample set that can be specified and be different for each
value of Gamma, when (M*H) files are allocated in the Helpers by means of
a DGPA Learning Algorithm. The script returns also the CDF of the
Average Network Weighted Delay, and the Average INITIALISATION and GAME 
complexity calculated in temrs of number of iterations.
           

----------------------------- DEPENDENCIES --------------------------------

INITIALIZE_Game_Of_DGPA (...)
PLAY_Game_Of_DGPA(... )
Network_A_W_D(...)

-------------------------------- INPUT  ----------------------------------
           
Results_Directory
           The Directory where figures are saved.
Reward_Type
           See INITIALIZE_Game and PLAY_Game functions.      
Network_Delays
           Matrix NxH+1 of delays between users and providers
Gamma_ZipF
           The range of parameters to try for the Zip-f probability ditrib.
SST_L
           Number of sampels for each leanring to infer some stats.
M
           Number of file memories in each helper.
H
           Number of helpers in the cell
N
           Number of users in the cell
S
           Library of files to cache
Ini_Number
           See NITIALIZE_Game_Of_DGPA
Resolution
           See PLAY_Game_Of_DGPA
Conv_Prob_Th
           See PLAY_Game_Of_DGPA
            
-------------------------------- OUTPUT  ----------------------------------

Failures
           

Avg_Net_AWD_DGPA
           Average Network Weighted Delay resulting from the Learned
           Caching

Iterations_DGPA_INIT
           Number of iterations needed to initialize the GAME
           
Iterations_DGPA_GAME
           Number of iterations needed to complete the GAME

Avg_Iter_GAME_DGPA
           Like Iterations_DGPA_GAME but avergaed across the statitcal
           sample SST_L

Avg_Iter_INIT_DGPA
           Like Iterations_DGPA_INIT but avergaed across the statitcal
           sample SST_L

% -----------------------------   CODE   ---------------------------------
%}

%Pre-allocate output variables for speed
Failures = zeros(N_Gammas,SST_L);
Avg_Net_AWD_DGPA = zeros(1,N_Gammas);
Avg_Iter_GAME_DGPA = zeros(1,N_Gammas);
Avg_Iter_INIT_DGPA = zeros(1,N_Gammas);
Iterations_DGPA_INIT = zeros(N_Gammas,SST_L);
Iterations_DGPA_GAME = zeros(N_Gammas,SST_L);


Net_AWD_DGPA = zeros(N_Gammas,SST_L);
N_Gammas = size(Gamma_ZipF,2);
Markers = {'--ob', '+m', '*c', 'xr', 'sg', 'db', '>b', 'pk', '--k', '--hr'};


% Cumulative Distribution Function of Network Average Weighted Delay
Figure_Handle_1 = figure;
Title_1 = ['CDF of NAWD for DGPA Placements: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_1 = fullfile(Results_Directory, ['CDF_NAWD_DGPA_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title_1);
xlabel('AWDs');
ylabel('P( NAWD <= x)');
grid on;
hold on;

time = 0;
Counter = 1;
Allocations_Left = N_Gammas*SST_L;

for g = 1:1:N_Gammas
    [ Popularities ] = Files_Popularities( S, Gamma_ZipF(g) );
    for i = 1:1:SST_L
        tstart=tic;
        
        try
            
           % Initialize the game
           [ Learning, Iterations_DGPA_INIT(g,i) ] = INITIALIZE_Game_Of_DGPA ( Network_Delays, ...
                                          M, Popularities, Reward_Type, Ini_Number );
           % Play Game                               
           [ Alloc_DGPA, Iterations_DGPA_GAME(g,i)] = PLAY_Game_Of_DGPA ( Network_Delays, ...
                           Popularities, Learning, Reward_Type, Resolution, Conv_Prob_Th );
                       
           % Calculate network average weighted delay for particular gamma
           % and this iteration out of SST_L iterations:
           Net_AWD_DGPA(g,i) =  Network_A_W_D( Network_Delays, S, Popularities, Alloc_DGPA );
           
        catch
            Failures(g,i) = 1;
        end
        
        time = time + toc(tstart);
        Avg_Time = time/Counter;
        Time_Left = Avg_Time*Allocations_Left;
        
        disp(['Learning. Time left: ' num2str(floor(Time_Left/60)), ...
            ' minutes and ' num2str(rem(Time_Left,60)) ' seconds.' ]);
        
        Counter = Counter + 1;
        Allocations_Left = Allocations_Left - 1;
    end
    
    Avg_Net_AWD_DGPA(g) = sum ( Net_AWD_DGPA(g,:) ) / SST_L; 
    
    % Average number of iterations needed to initialize and to play
    Avg_Iter_GAME_DGPA(g) = sum ( Iterations_DGPA_GAME(g,:) ) / SST_L; 
    Avg_Iter_INIT_DGPA(g) = sum ( Iterations_DGPA_INIT(g,:)) / SST_L; 
    
    % Empirical (Kaplan-Meier) cumulative distribution function.
    [CDF_AWD_DGPA, Set_of_AWD_DGPA] = ecdf( Net_AWD_DGPA(g,:) );
    
    % Plot CDF
    plot(Set_of_AWD_DGPA, CDF_AWD_DGPA, Markers{g}); 
    hold on;
end

legend(['Gamma 0.1','Gamma 0.2','Gamma 0.3','Gamma 0.4','Gamma 0.5',...
        'Gamma 0.6','Gamma 0.7','Gamma 0.8','Gamma 0.9','Gamma 1']);
    
cd(Results_Directory);
saveas(Figure_Handle_1,Fig_File_Name_1);


% Average Network AWD
Figure_Handle_2 = figure;
plot(Gamma_ZipF, Avg_Net_AWD_DGPA,'k'); 
Title_2 = ['Average NAWD for DGPA Placements: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_2 = fullfile(Results_Directory, ['ANAWD_DGPA_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title_2);
xlabel('Gamma');
ylabel('Average NAWD');
grid on;

cd(Results_Directory);
saveas(Figure_Handle_2,Fig_File_Name_2);

% AVERAGE COMPLEXITY INITIALISATION
Figure_Handle_3 = figure;
plot(Gamma_ZipF, Avg_Iter_INIT_DGPA,'--xr'); 
hold on;
plot(Gamma_ZipF, Avg_Iter_GAME_DGPA, 'ob'); 
legend('INIT_DGPA','GAME_DGPA');
Title_3 = ['Average Complexity ' num2str(SST_L) ' DGPA Placements: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_3 = fullfile(Results_Directory, ['Avg_Complexity_DGPA_' num2str(SST_L) '_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title_3);
xlabel('Gamma');
ylabel('Iterations');
grid on;
hold on;

cd(Results_Directory);
saveas(Figure_Handle_3,Fig_File_Name_3);

end

