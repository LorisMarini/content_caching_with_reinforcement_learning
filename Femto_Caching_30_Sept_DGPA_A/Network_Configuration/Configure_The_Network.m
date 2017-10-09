function [ Network_Parameters ] = Configure_The_Network( N_Networks, Diameter, Radius_Protected_Area, Step_Size, H, N, Alpha)
% ATTENTION: 
% Make sure that you run this function only once as it may overwrite the
% existing network configuration and cause a mess. Here it is an example 
% of how to call it to generate 100 different Networks.
% Configure_The_Network( 100, Diameter,Radius_Protected_Area,Step_Size, H, N, Alpha);  

Network_Config_Direcotry = 'C:\Users\lmar1564\Documents\MATLAB\FemtoChaching\Network_Configuration\Configurations';

Network_Parameters = struct('Number_Of_Netoworks_Generated',N_Networks, 'Diameter',Diameter, ...
                     'Radius_Protected_Area',Radius_Protected_Area, 'Step_Size',Step_Size, ...
                     'N_Of_Helpers',H, 'N_Of_Users', N, 'Propagation_Loss_Alpha', Alpha );         
Parameters_File_Name = ['Network_Parameters_' num2str(N_Networks) '_Random_Networks_H' num2str(H) '_N' num2str(N) '.mat'];
save( fullfile( Network_Config_Direcotry ,Parameters_File_Name),  'Network_Parameters');


for n = 1:1:N_Networks
    
    [ Distances, Fig_Handle ] = Place_Users( Diameter, Radius_Protected_Area, Step_Size, H,N );
    Network_Delays = Network_Normalised_Delays( Distances, Alpha );
    
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
end

close all;
end

