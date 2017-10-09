% THIS IS A STABILITY ANALYSIS.
% WE ANALYSE THE STATISTICS OF THE FIRST ORDER, RANGE AND EXPECTATION.
% WE WOULD LIKE TO HAVE THE EXPECTATION AS CLOSE AS POSSIBLE TO THE
% Benchmark (GREEDY) and THE RANGE (DeltaRelative) AS SMALL AS POSSIBLE.

% Statistical Test for Game_DGPA_FC_Function
% Author: Loris Marini
% Version: 1.0.1
% Date: 18/09/2014

clear all;
close all;
format shortE;

% FEMTO CACHING TEST PARAMETERS

Diameter = 1000;
Radius_Protected_Area = 15;
Step_Size = 0.1;
H = 5;
N = 10;
Alpha = 3;
Gamma_ZipF = 0.4;
M = 3;                                  % Caching capability of each Helper.
S = 1:1:H*M;                            % Set of files we can off-load.
Resolution = 1;                         % DGPA Learning Resolution
Conv_Prob_Th = 0.9;                     % DGPA Convergence Threshold
Ini_Number = 10;                        % DGPA Initialisation Minimum Selection
SST = 5;

Rewards = { 'Best_File_Based_Reward', 'Weighted_Delay_Based_Reward','Average_Weighted_Delay_Based_Reward'};
Reward_Type = Rewards{2};

Configure_The_Network( 1, Diameter,Radius_Protected_Area,Step_Size, H, N, Alpha);
load('Network_Delays_1');

[ Popularities ] = Files_Popularities( S, Gamma_ZipF );

WD_DGPA_Samples = zeros(SST,N);         % SST samples of the Weighted Delays computed with DGPA Learning Automata
A_WD_DGPA_Samples = zeros(SST,N);       % SST samples of the Average Weighted Delays computed with DGPA Learning Automata
A_WD_Greedy_Samples = zeros(SST,N);     % SST samples of the Avergae Weighted Delays computed with Greedy Algorithm.
A_WD_Greedy_Fast_Samples = zeros(SST,N);% SST samples of the Avergae Weighted Delays computed with Fast Greedy Algorithm.

for i=1:1:SST

    [ Learning ] = INITIALIZE_Game_Of_DGPA( Network_Delays, M, Popularities, Reward_Type, Ini_Number );
    [ Learning_Fast] = Fast_Game_DGPA_INITIALISATION ( Network_Delays, M, Popularities, Reward_Type, Ini_Number );

    [ ~, ~,~, Conv_Delays_DGPA, ~] ...
        = PLAY_Game_Of_DGPA( Network_Delays, Popularities, Learning, Reward_Type, Resolution, Conv_Prob_Th );
        
    [ ~, ~, ~,Conv_Delays_DGPA_FAST, ~] ...
        = PLAY_Game_Of_DGPA( Network_Delays, Popularities, Learning_Fast, Reward_Type, Resolution, Conv_Prob_Th );
    
    A_WD_DGPA_Samples(i,:) = Conv_Delays_DGPA;
    A_WD_DGPA_Samples_FAST(i,:) = Conv_Delays_DGPA_FAST;
%{
    [ OA_Greedy_Fast, A_WD_Greedy_Fast ] = Greedy_Placement_Fast( Network_Delays, M, Popularities ); 
    [ OA_Greedy, A_WD_Greedy ] = Greedy_Placement( Network_Delays, M, Popularities );
    A_WD_Greedy_Samples(i,:) = A_WD_Greedy;
    A_WD_Greedy_Fast_Samples(i,:) = A_WD_Greedy_Fast;
    %}
end


Range_DGPA = abs( max(A_WD_DGPA_Samples) - min(A_WD_DGPA_Samples));
Expectation_DGPA = sum(A_WD_DGPA_Samples,1)./SST;
R_DGPA = (Range_DGPA./Expectation_DGPA)*100;

Range_DGPA_FAST = abs( max(A_WD_DGPA_Samples_FAST) - min(A_WD_DGPA_Samples_FAST));
Expectation_DGPA_FAST = sum(A_WD_DGPA_Samples_FAST,1)./SST;
R_DGPA_FAST = (Range_DGPA_FAST./Expectation_DGPA_FAST)*100;

%{
Range_Greedy = abs( max(A_WD_Greedy_Samples) - min(A_WD_Greedy_Samples));
Expectation_Greedy = sum(A_WD_Greedy_Samples,1)./SST;
R_Greedy = (Range_Greedy./Expectation_Greedy)*100;

Range_Greedy_Fast = abs( max(A_WD_Greedy_Fast_Samples) - min(A_WD_Greedy_Fast_Samples));
Expectation_Greedy_Fast = sum(A_WD_Greedy_Fast_Samples,1)./SST;
R_Greedy_Fast = (Range_Greedy_Fast./Expectation_Greedy_Fast)*100;
%}
stop=1;
