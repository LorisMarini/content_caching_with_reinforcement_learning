% -------------------------------------------------------------------------
% This Script has been designed to test the performance of the proposed
% Decentralised Learning scheme based on DGPA Learning Automata, against the 
% existing Greedy Centraliased Algorithm. 

clear all; 
close all;

%% Directory Set-Up
Res_Dir ='C:\Users\lmar1564\Dropbox\Research USYD\Femto_Caching_Simulations\FemtoChaching_4_Oct_Cooperative_Learners\Large Scale Performance\Results';
Network_Config_Direcotry = 'C:\Users\lmar1564\Dropbox\Research USYD\Femto_Caching_Simulations\FemtoChaching_4_Oct_Cooperative_Learners\Large Scale Performance\Networks';
Final_Res_Dir = 'C:\Users\lmar1564\Dropbox\Research USYD\Femto_Caching_Simulations\FemtoChaching_4_Oct_Cooperative_Learners\Large Scale Performance\Final Results';

%% Network Parameters
Diameter = 1000;                  % Length of the cell edge (square model)
Radius_Protected_Area = 15;       % Area close to MBS and Helpers where users cannot be placed.
Step_Size = 0.1;                  % Unit of distance along x-axis or y-axis (meters)
H = 12;                           % Number of Helpers
N = 100;                          % Number of Users
Alpha = 3;                        % Exponential Path-Loss Coefficient
M = 3;                            % Caching capability of each Helper;
Gamma_ZipF = [0.1, 0.5, 1];       % Discounted Rate (Zip-F Popularity Distribution);
N_Gammas = size(Gamma_ZipF,2);   

%% Learning Parameters
Resolution = 1;                   % DGPA Learning Resolution;
Conv_Prob_Th = 0.999;             % DGPA Convergence Threshold;
Ini_Number = 5;                   % DGPA Initialisation Minimum Selection;
Reward_Type = 'Average_Weighted_Delay_Based_Reward';             
Learning_Setup = 'Cooperative_Shuffle';

%% Parameters for Statistical Analysis 

N_Networks = 1000;                % Number of Networks (Different User Placement, same Helpers Location)
SST = 1;                          % Number of symulations per gamma, per network

%% Variables Pre-Allocation

Avg_NAWD_Greedy = zeros(N_Networks,N_Gammas);
Avg_Net_AWD_DGPA_CS = zeros(N_Networks,N_Gammas);
Avg_Iter_GAME_DGPA_CS = zeros(N_Networks,N_Gammas);
Avg_Iter_INIT_DGPA_CS = zeros(N_Networks,N_Gammas);
NAWD_Uncached = zeros(N_Networks,N_Gammas);
Tot_Failures = zeros(N_Gammas, SST);

%% Performance Analysis
for n = 1:1:N_Networks
    
    Net_Number = n;
    
    [ Distances, Fig_Handle ] = Place_Users( Diameter, Radius_Protected_Area, Step_Size, H, N );
    Network_Delays = Network_Normalised_Delays( Distances, Alpha );
    
    %{
    % SAVE DISTANCES
    Matrix_Distances_Name = ['Distances_' num2str(n) '_of_' num2str(N_Networks) '_with_H' num2str(H) '_N' num2str(N) '.mat'];
    Distances_File_Name = fullfile( Network_Config_Direcotry ,Matrix_Distances_Name);
    save( Distances_File_Name, 'Distances');
    % SAVE DELAYS
    Matrix_Delays_Name = ['Network_Delays_' num2str(n) '_of_' num2str(N_Networks) '_with_H' num2str(H) '_N' num2str(N) '.mat'];
    Delays_File_Name = fullfile( Network_Config_Direcotry ,Matrix_Delays_Name);
    save( Delays_File_Name, 'Network_Delays');
    % SAVE PICTURE
    Figure_Name = ['Network_Scenario_' num2str(n) '_of_' num2str(N_Networks) '_with_H' num2str(H) '_N' num2str(N) '.fig'];
    Figure_File_Name = fullfile(Network_Config_Direcotry, Figure_Name);
    figure(Fig_Handle);
    saveas(Fig_Handle, Figure_File_Name);
    %}
    H = size(Network_Delays, 2) - 1;          % Number of Helpers in the cell;
    N = size(Network_Delays, 1);              % Number of Users in the cell;
    F = H*M;                                  % Max files we can off-load;
    S = 1:1:F;                                % Set of files we can off-load;
    try
        [ NAWD_Uncached(n,:) ] = Large_Scale_Performance_Uncached( Net_Number, N_Networks, Res_Dir,Network_Delays, Gamma_ZipF, M, H, N, S  );
        close all
    catch
        NAWD_Uncached(n,:) = NAWD_Uncached(n-1,:);
        UC_Errors = UC_Errors +1;
    end
    try
        [ Avg_NAWD_Greedy(n,:), ~ ]= Large_Scale_Performance_Greedy_Placement(Net_Number, N_Networks, Res_Dir, Network_Delays, Gamma_ZipF, SST, M, H, N, S  );
        close all
    catch
        Avg_NAWD_Greedy(n,:) = Avg_NAWD_Greedy(n-1,:);
        G_Errors = G_Errors +1;
    end
    try
        [ Failures_CS, Avg_Net_AWD_DGPA_CS(n,:), Avg_Iter_GAME_DGPA_CS(n,:), Avg_Iter_INIT_DGPA_CS(n,:),~, ~ ] = Large_Scale_Performance_DGPA(Net_Number,N_Networks, Learning_Setup, Res_Dir, Reward_Type, Network_Delays, Gamma_ZipF, SST, M, H, N, S, Ini_Number, Resolution, Conv_Prob_Th );
        Tot_Failures = Tot_Failures + Failures_CS;
        close all
    catch
        Avg_Net_AWD_DGPA_CS(n,:) = Avg_Net_AWD_DGPA_CS(n-1,:);
        L_Errors = L_Errors +1;
    end
end

Avg_NAWD_Uncached = sum(NAWD_Uncached,1)./size(NAWD_Uncached,1);
Avg_Avg_NAWD_Greedy = sum(Avg_NAWD_Greedy,1)./size(Avg_NAWD_Greedy,1);
Avg_Avg_Net_AWD_DGPA_CS = sum(Avg_Net_AWD_DGPA_CS,1)./size(Avg_Net_AWD_DGPA_CS,1);
Avg_Avg_Iter_GAME_DGPA_CS = sum(Avg_Iter_GAME_DGPA_CS,1)./size(Avg_Iter_GAME_DGPA_CS,1);
Avg_Avg_Iter_INIT_DGPA_CS = sum(Avg_Iter_INIT_DGPA_CS,1)./size(Avg_Iter_INIT_DGPA_CS,1);

Figure_Handle_1 = figure;
hold on;
bar([Avg_Avg_NAWD_Greedy',Avg_Avg_Net_AWD_DGPA_CS', Avg_NAWD_Uncached']);
legend('Greedy','DGPA');
Title_1 = ['ANWD Comparison with ' num2str(N_Networks) ' Networks and ' num2str(SST) ' samples: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_1 = fullfile(Final_Res_Dir, ['ANWD_Comparison_'  num2str(N_Networks) '_Networks_' num2str(SST) '_Samples_H_' num2str(H) '_M_' num2str(M) '_N_' num2str(N)  '.fig']);
title(Title_1);
xlabel('Gamma');
ylabel('Average ANWD');
grid on;
saveas(Figure_Handle_1,Fig_File_Name_1);

%%
Markers = {'--ob', '+m', '*c', 'xr', 'sg', 'db', '>b', 'pk', '--k', '--hr'};
Figure_Handle_2 = figure;
hold on;
Leg =[];
for g=1:1:N_Gammas
    [CDF_Uncached, Sets_Uncached] = ecdf( NAWD_Uncached(:,g) );
    [CDF_Learning, Sets_Learning] = ecdf( Avg_Net_AWD_DGPA_CS(:,g) );
    [CDF_Greedy, Sets_Greedy] = ecdf( Avg_NAWD_Greedy(:,g) );
    plot(Sets_Uncached, CDF_Uncached, Markers{((g-1)*3 +1)});
    plot(Sets_Learning, CDF_Learning, Markers{((g-1)*3 +2)});
    plot(Sets_Greedy, CDF_Greedy, Markers{((g-1)*3 +3)});
    lgnd{((g-1)*3 +1)} = ['Uncached ' num2str(Gamma_ZipF(g))];
    lgnd{((g-1)*3 +2)} = ['Learning ' num2str(Gamma_ZipF(g))];
    lgnd{((g-1)*3 +3)} = ['Greedy ' num2str(Gamma_ZipF(g))]; 
end 
Title_2 = ['CDF Comparison with ' num2str(N_Networks) ' Networks and ' num2str(SST) ' samples: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_2 = fullfile(Final_Res_Dir, ['CDF_Comparison_'  num2str(N_Networks) '_Networks_' num2str(SST) 'samples_H_' num2str(H) '_M_' num2str(M) '_N_' num2str(N)  '.fig']);
title(Title_2);
xlabel('NAWD');
ylabel('CDF of NAWD');
grid on;
saveas(Figure_Handle_2,Fig_File_Name_2);
legend(lgnd);
save(['Large_Scale_Simulatios_' num2str(N_Networks) '_Networks_' num2str(SST) '_Samples_H_' num2str(H) '_M_' num2str(M) '_N_' num2str(N) '.mat']);


