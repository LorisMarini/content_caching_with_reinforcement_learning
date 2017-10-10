function [ Avg_Net_AWD_Greedy, Net_AWD_Greedy ] = Performance_Greedy_Placement( Results_Directory, Network_Delays, Gamma_ZipF, SST_G, M, H, N, S  )

% This Function sets out to calculate the Average of the Average Weighted
% Delay of the users for a set of Gammas (Zip-f distribution parameter) and
% a statistical sample set that can be specified and be different for each
% value of Gamma, when (M*H) files are allocated in the Helpers by means of
% a Greedy placement Algorithm. The script returns also the CDF of the
% Average Network Weighted Delay, and the complexity calculated in temrs of
% number of iterations.

% Author: LORIS MARINI
% Date: 26/09/2014
% Version: 1.0

N_Iterations = H*nchoosek(M*H,M);
N_Gammas = size(Gamma_ZipF,2);
Net_AWD_Greedy = zeros(N_Gammas,SST_G);
Avg_Net_AWD_Greedy = zeros(1,N_Gammas);
Markers = {'--ob', '+m', '*c', 'xr', 'sg', 'db', '>b', 'pk', '--k', '--hr'};

% CDFs Network Average Weighted Delay
Figure_Handle_3 = figure;
Title_3 = ['CDF of NAWD with ' num2str(SST_G) ' Greedy Placements: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_3 = fullfile(Results_Directory, ['CDF_NAWD_Greedy_' num2str(SST_G) '_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title_3);
xlabel('AWDs');
ylabel('P( NAWD <= x)');
grid on;
hold on;

T = 0;
Counter = 1;
Allocations_Left = N_Gammas*SST_G;

for g = 1:1:N_Gammas

    [ Popularities ] = Files_Popularities( S, Gamma_ZipF(g) );
    
    for i = 1:1:SST_G
        tstart = tic;
        [ ALL_Greedy, ~ ] = Greedy_Placement_Fast( Network_Delays, M, Popularities );
        Net_AWD_Greedy(g,i) =  Network_A_W_D( Network_Delays, S, Popularities, ALL_Greedy );
        T = T + toc(tstart);
        Avg_Time = T/Counter;
        Time_Left = Avg_Time*Allocations_Left;
        disp(['Greedy Placement. Time left: ', num2str(floor(Time_Left/60)), ' minutes and ' num2str(rem(Time_Left,60)) ' seconds.' ]);
        Counter = Counter + 1;
        Allocations_Left = Allocations_Left - 1;
    end
    Avg_Net_AWD_Greedy(g) = sum ( Net_AWD_Greedy(g,:) ) / SST_G; 
    [CDF_AWD_Greedy, Set_of_AWD_Greedy] = ecdf( Net_AWD_Greedy(g,:) );
    plot(Set_of_AWD_Greedy,CDF_AWD_Greedy, Markers{g}); 
    hold on;
end

legend('Gamma 0.1','Gamma 0.2','Gamma 0.3','Gamma 0.4','Gamma 0.5','Gamma 0.6','Gamma 0.7','Gamma 0.8','Gamma 0.9','Gamma 1');
saveas(Figure_Handle_3,Fig_File_Name_3);

% Average Network AWD
figure;
plot(Gamma_ZipF, Avg_Net_AWD_Greedy,'k'); 
Title_4 = ['Average NAWD with ' num2str(SST_G) ' Greedy Placements: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_4 = fullfile(Results_Directory, ['Avg_NAWD_Greedy_' num2str(SST_G) '_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title_4);
xlabel('Gamma');
ylabel('Average NAWD');
grid on;
saveas(gcf,Fig_File_Name_4);

% COMPLEXITY
figure;
plot(Gamma_ZipF, N_Iterations*ones(1,N_Gammas)); 
Title_5 = ['Average Complexity ' num2str(SST_G) ' Greedy Placements: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_5 = fullfile(Results_Directory, ['Avg_Complexity_Greedy_' num2str(SST_G) '_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title_5);
xlabel('Gamma');
ylabel('Average NAWD');
grid on;
saveas(gcf,Fig_File_Name_5);
end



