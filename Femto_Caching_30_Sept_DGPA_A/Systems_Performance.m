% SYSTEM PERFORMANCE
% We set out to calculate the average weighted delay among all the users
% in the network, for different types of algorithms and different
% popularity distributions. 
clear all;
close all;

Rewards = { 'Best_File_Based_Reward', 'Weighted_Delay_Based_Reward',...           
            'Average_Weighted_Delay_Based_Reward'};

% LOAD THE NETWORK
load('Network_Delays_1_of_10_with_H4_N100');     % We focus on the first configuration;
H = size(Network_Delays, 2) - 1;                  % Number of Helpers in the cell;
N = size(Network_Delays, 1);                      % Number of Users in the cell;

M = 3;                                            % Caching capability of each Helper;
F = H*M;                                          % Max files we can off-load;
S = 1:1:F;                                        % Set of files we can off-load;
Gamma_ZipF = 0.1 : 0.1 : 1;                       % Different Popularity Distributions;
N_Gammas = size(Gamma_ZipF,2);

Resolution = 1;                                   % DGPA Learning Resolution;
Conv_Prob_Th = 0.999;                             % DGPA Convergence Threshold;
Ini_Number = 10;                                  % DGPA Initialisation Minimum Selection;
Reward_Number = 2;                                % Reward Strategy of the DGPA
Reward_Type = Rewards{Reward_Number};             

SST_Rnd = 1000;            % Number of points among which we average for the Random Case;
SST_G = 100;               % Number of points among which we average for the Greedy Case;
SST_L = 10;               % Number of points among which we average for the Learning Case;

Results_Directory ='C:\Users\lmar1564\Documents\MATLAB\FemtoChaching\Performance\Plots';

%{
[ Net_AWD_Uncached ] = Performance_Uncached(Results_Directory, Network_Delays, Gamma_ZipF, M, H, N, S  );
[ Avg_Net_AWD_Greedy, Net_AWD_Greedy ]= Performance_Greedy_Placement(Results_Directory, Network_Delays, Gamma_ZipF, SST_G, M, H, N, S  );
[ Avg_Net_AWD_Random, Net_AWD_Random ] = Performance_Random_Placement(Results_Directory, Network_Delays, Gamma_ZipF, SST_Rnd, M, H, N, S  );
%}

[ Failures, Avg_Net_AWD_DGPA, Avg_Iter_GAME_DGPA, Avg_Iter_INIT_DGPA, Iterations_DGPA_INIT, Iterations_DGPA_GAME ] = ...
    Performance_Learning(Results_Directory, Reward_Type, Network_Delays, Gamma_ZipF, SST_L, M, H, N, S, Ini_Number, Resolution, Conv_Prob_Th );

%{
[ Net_AWD_Best_Random ] = Performance_Best_Random_Placement(Results_Directory, Network_Delays, Gamma_ZipF, 5*ones(1,10), M, H, N, S  );
%}