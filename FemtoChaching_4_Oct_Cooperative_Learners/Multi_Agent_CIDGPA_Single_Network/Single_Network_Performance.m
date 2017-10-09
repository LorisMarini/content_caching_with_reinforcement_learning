% ---------------- FEMTO CACHING: SYSTEM PERFORMANCE-----------------------
% This script sets out to compare the Average Weighted Delay among the users 
% of a cell in a cellular network, for a set of Gammas (Zip-f distribution 
% parameter) when (M*H) files are allocated in the Helpers by means of
% different algorithms. We compare the Uncached case with the Random,
% Greedy, Learning and Best Random approach(*). 

% Author: LORIS MARINI
% Date: 20/10/2014
% Version: 1.0.2

%(*) Not included in the ICC2015 Publication.

clear all; 
close all;
Results_Directory ='E:\DATA';

load('Network_Delays_1_of_10_with_H12_N100');     % We focus on the first configuration;
H = size(Network_Delays, 2) - 1;          % Number of Helpers in the cell;
N = size(Network_Delays, 1);              % Number of Users in the cell;
M = 3;                                    % Caching capability of each Helper;
F = H*M;                                  % Max files we can off-load;
S = 1:1:F;                                % Set of files we can off-load;
Gamma_ZipF = 0.1 : 0.1 : 1;               % Different Popularity Distributions;
N_Gammas = size(Gamma_ZipF,2);

%% Learning Algorithm Parameters

Rewards = { 'Best_File_Based_Reward', ...
            'Weighted_Delay_Based_Reward', ...
            'Average_Weighted_Delay_Based_Reward'};
        
Possible_Setups = {'Single',...
                   'Cooperative_Hard', ...
                   'Cooperative_Shuffle'};
               
Resolution = 1;                           % Learning Resolution;
Conv_Prob_Th = 0.999;                     % Convergence Threshold;
Ini_Number = 5;                           % Initialisation Minimum Selection;
Learning_Setup = Possible_Setups{3};      % The way Learners Search the file.
Reward_Type = Rewards{3};                 % Reward Function for the game of CI-DGPA

%% Statistical Analysis Parameters

SST_Rnd = 1000;   % Number of points among which we average for the Random Case;
SST_G = 30;       % Number of points among which we average for the Greedy Case;
SST_L = 30;       % Number of points among which we average for the Learning Case;

%% Single Network Performance (Uncomment To Run)

% UNCACHED
%[ Net_AWD_Uncached ] = Performance_Uncached(Results_Directory,Network_Delays, Gamma_ZipF, M, H, N, S  );

% RANDOM
%[ Avg_Net_AWD_Random, Net_AWD_Random ] = Performance_Random_Placement(Results_Directory, Network_Delays, Gamma_ZipF, SST_Rnd, M, H, N, S  );

% GREEDY
%[ Avg_Net_AWD_Greedy, Net_AWD_Greedy ]= Performance_Greedy_Placement(Results_Directory, Network_Delays, Gamma_ZipF, SST_G, M, H, N, S  );

% LEARNING 1
%[ Failures_S, Avg_Net_AWD_DGPA_S, Avg_Iter_GAME_DGPA_S, Avg_Iter_INIT_DGPA_S,Iterations_DGPA_INIT_S, Iterations_DGPA_GAME_S ] = Performance_Learning(...
%    Possible_Setups{1}, Results_Directory, Reward_Type, Network_Delays, Gamma_ZipF, SST_L, M, H, N, S, Ini_Number, Resolution, Conv_Prob_Th );

% LEARNING 1
%[ Failures_CH, Avg_Net_AWD_DGPA_CH, Avg_Iter_GAME_DGPA_CH, Avg_Iter_INIT_DGPA_CH,Iterations_DGPA_INIT_CH, Iterations_DGPA_GAME_CH ] = Performance_Learning(...
%    Possible_Setups{2}, Results_Directory, Reward_Type, Network_Delays, Gamma_ZipF, SST_L, M, H, N, S, Ini_Number, Resolution, Conv_Prob_Th );

% LEARNING 1
[ Failures_CS, Avg_Net_AWD_DGPA_CS, Avg_Iter_GAME_DGPA_CS, Avg_Iter_INIT_DGPA_CS,Iterations_DGPA_INIT_CS, Iterations_DGPA_GAME_CS ] = Performance_Learning(...
    Possible_Setups{3}, Results_Directory, Reward_Type, Network_Delays, Gamma_ZipF, SST_L, M, H, N, S, Ini_Number, Resolution, Conv_Prob_Th );

% BEST RANDOM                          
% [ Net_AWD_Best_Random ] = Performance_Best_Random_Placement(Results_Directory,Network_Delays, Gamma_ZipF, Avg_Iter_GAME_DGPA, M, H, N, S  );
      

% save('Best_Learning_LCS_Nini5_RW3_1_of_10_with_H12_N100.mat');


