function [ Net_AWD_Uncached ] = Performance_Uncached( Results_Directory, Network_Delays, Gamma_ZipF, M, H, N, S  )

N_Gammas = size(Gamma_ZipF,2);
Length_Popularities = M*H;
NetAWD_Uncached = zeros(1,N_Gammas);

for g = 1:1:N_Gammas
    [ Popularities ] = Files_Popularities( S, Gamma_ZipF(g) );
    Net_AWD_Uncached(g) =  sum( diag( repmat(Popularities,N,1)* ( repmat(Network_Delays(:,end),1,Length_Popularities))'  )) /N;
end

% Average Network AWD
figure;
plot(Gamma_ZipF, Net_AWD_Uncached,'k'); 
Title = ['NAWD No Caching: H=' num2str(H)  ' M=' num2str(M) ' N=' num2str(N)];
Fig_File_Name = fullfile(Results_Directory, ['NAWD_Uncached_N_' num2str(N) '_H' num2str(H)  '_M' num2str(M) '.fig']);
title(Title);
xlabel('Gamma');
ylabel('Average NAWD');
grid on;
saveas(gcf,Fig_File_Name);
end

