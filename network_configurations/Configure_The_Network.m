function [ Network_Parameters ] = Configure_The_Network( N_Networks, Diameter, Radius_Protected_Area, Step_Size, H, N, Alpha)

%{
-------------------------   AUTHORSHIP  -------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

-------------------------   DESCRIPTION   -------------------------

This function creates N_Networks different instances of a radio cell. If the 
function is run more than once there are two possibilities:

1. The Network_Parameters differs from those defined previously. The function 
   creates a subdirectory with a random name and saves the data there.

2. The Network_Parameters is identical to the one already defined. The
   function does not do anyhting.
 

------------------------- INPUT PARAMETERS -------------------------

-- N_Networks --
Number of different networks to simulate

-- Diameter -- 
Diameter of the cell in meters.

-- Radius_Protected_Area --
Radius of clear space around each user (Necessary to avoid ...)

-- Step_Size --
Distance between two points on the cell. Must be <= 1. The ratio Diameter/Step_Size 
is equal to the number of points along the diameter of the cell.

-- H -- 
Number of Helpers in the cell.

-- N --
Number of users in the cell.

-- Alpha --
The free space absorption coefficient 

------------------------- OUTPUT PARAMETERS -------------------------

-- Parameters -- 
A struct containing all the relevant parameters for the function call.

------------------------- EXAMPLE OF CALL -----------------------

% Number of radio cells to deploy
Nrc = 2;

% Cell diameter [m]
Diam = 1000;

% Clear space around each helper [m]:
CS = 10;

% Step size
SS = 0.1

% Number of helpers in the network
Nh = 20;

% Number of users in the network
Nu = 100;

% A propagation Loss of 1 [km^-1]
Ap = 1;

[ Parameters ] = Configure_The_Network( Nrc, Diam, CS, SS, Nh, Nu, Ap)

%}
% ----------------------------   CODE     --------------------------

% Add the communications toolbox to path. This is a known issue in Matlab
% 2017a. See this post for more info: 

% https://au.mathworks.com/matlabcentral/answers/350491-how-to-set-the-correct-path-
% for-communications-toolbox-for-rtl-sdr-support-package-at-startup
        
addpath(fullfile(matlabroot, 'toolbox', 'comm', 'comm'), '-end')

% Define the path to the directory where you store the network configurations:
% path_to_networks = 'C:\Users\lmar1564\Documents\MATLAB\FemtoChaching\Network_Configuration\Configurations';
path_to_networks = '/home/loris/Desktop/test';

                        
% Create a struct with the name and value of all relevant parameters used
% in this function call:
Network_Parameters = struct('Number_Of_Netoworks_Generated',N_Networks, 'Diameter',Diameter, ...
                     'Radius_Protected_Area',Radius_Protected_Area, 'Step_Size',Step_Size, ...
                     'N_Of_Helpers',H, 'N_Of_Users', N, 'Propagation_Loss_Alpha', Alpha ); 
            
                 
% Define the name with which you save the Network_Parameters:             
parameters_saveas = ['Network_Parameters_' num2str(N_Networks),...
                        '_Random_Networks_H' num2str(H) '_N' num2str(N) '.mat'];
                
                    
cd(path_to_networks);            
search = dir(parameters_saveas);

if ~isempty( search )
    % Check they are identical
    l = load(parameters_saveas);
    Loaded_param = l.Network_Parameters;
    
    if isequal(Loaded_param, Network_Parameters)
        % File exists and has identical parameters. Do nothing
        disp('Configuration already exists.');
        return;
        
    elseif Loaded_param ~= Network_Parameters
        % perform the full simulation but save data in subdirectory that
        % starts with the idenfifiers of parameters_saveas plus a random
        % string of integers
        
        random_string = num2str(randi(9,1,7));
        random_string(random_string==' ') = [];
        
        subdir_name = ['N_Net_' num2str(N_Networks),'H' num2str(H) '_N' num2str(N),'_id_',random_string];
        % ????        
        
    end
else
    % Do nothing here and continue execution.
end

            
% Save to disk the network parameters in the 'path_to_networks' with the 'parameters_saveas' name:                    
save( fullfile( path_to_networks ,parameters_saveas),  'Network_Parameters');


for n = 1:1:N_Networks
    
    % Call the function Place_Users to randomly place users in a cell
    [ distances, fh ] = Place_Users( Diameter, Radius_Protected_Area, Step_Size, H,N );
        
    % Save distances matrix to file 
    Matrix_Distances_Name = ['Distances_' num2str(n) '_of_' num2str(N_Networks) '_with_H' num2str(H) '_N' num2str(N) '.mat'];
    save( fullfile( path_to_networks ,Matrix_Distances_Name), 'distances');
    
    
    Network_Delays = Network_Normalised_Delays( distances, Alpha );
        
    % Save Matrix of delays
    Matrix_Delays_Name = ['Network_Delays_' num2str(n) '_of_' num2str(N_Networks) '_with_H' num2str(H) '_N' num2str(N) '.mat'];
    save( fullfile( path_to_networks ,Matrix_Delays_Name), 'Network_Delays');
    
    % Save the figure representing the cell with helpers and users
    Figure_Name = ['Network_Scenario_' num2str(n) '_of_' num2str(N_Networks) '_with_H' num2str(H) '_N' num2str(N) '.fig'];
    Figure_File_Name = fullfile(path_to_networks, Figure_Name);
    
    figure(fh);
    saveas(fh, Figure_File_Name);
end

close all;
end

