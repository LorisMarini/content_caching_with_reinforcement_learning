function [ Avg_Net_AWD_Random, Net_AWD_Random ] = Performance_Random_Placement( Results_Directory, Network_Delays, Gamma_ZipF, SST_Rnd, M, H, N, S  )

N_Gammas = size(Gamma_ZipF,2);
time = 0;
Counter = 1;
Allocations_Left = N_Gammas*SST_Rnd;
Net_AWD_Random = zeros(N_Gammas,SST_Rnd);

Markers = {'--ob', '+m', '*c', 'xr', 'sg', 'db', '>b', 'pk', '--k', '--hr'};

for g = 1:1:N_Gammas
    
    [ Popularities ] = Files_Popularities( S, Gamma_ZipF(g) );
    for i = 1:1:SST_Rnd
        tstart=tic;
        [ ALL_Random ] = Random_Placement( M,H );
        Net_AWD_Random(g,i) = Network_A_W_D( Network_Delays, S, Popularities, ALL_Random );  
        time = time + toc(tstart);
        Avg_Time = time/Counter;
        Time_Left = Avg_Time*Allocations_Left;
        disp(['Random Placement. Time left: ' num2str(floor(Time_Left/60)), ' minutes and ' num2str(rem(Time_Left,60)) ' seconds.' ]);
        Counter = Counter + 1;
        Allocations_Left = Allocations_Left - 1;
    end
    Avg_Net_AWD_Random(g) = sum ( Net_AWD_Random(g,:) ) / SST_Rnd;
    [CDF_AWD_Random(:,g), Set_of_AWD_Random(:,g)] = ecdf(Net_AWD_Random(g,:));
end

% Average Network AWD
Figure_Handle_1 = figure;
plot(Gamma_ZipF,Avg_Net_AWD_Random,'k'); 
Title = ['Average NAWD with ' num2str(SST_Rnd) ' Random Placements: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_1 = fullfile(Results_Directory, ['Avg_NAWD_Rnd_' num2str(SST_Rnd) '_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title);
xlabel('Gamma');
ylabel('Average NAWD');
grid on;
saveas(Figure_Handle_1,Fig_File_Name_1);

% CDFs Network Average Weighted Delay
Figure_Handle_2 = figure;
for g = 1:1:N_Gammas
    plot(Set_of_AWD_Random(:,g),CDF_AWD_Random(:,g), Markers{g});
    hold on;
end
Title_2 = ['CDF of NAWD with ' num2str(SST_Rnd) ' Random Placements: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name_2 = fullfile(Results_Directory, ['CDF_NAWD_Rnd_' num2str(SST_Rnd) '_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title_2);
xlabel('AWDs');
ylabel('P( NAWD <= x)');
grid on;
legend('Gamma 0.1','Gamma 0.2','Gamma 0.3','Gamma 0.4','Gamma 0.5','Gamma 0.6','Gamma 0.7','Gamma 0.8','Gamma 0.9','Gamma 1');
saveas(Figure_Handle_2,Fig_File_Name_2);

end



