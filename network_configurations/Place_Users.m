function [ distances, Fig_Handle ] = Place_Users( Diameter,Radius_Protected_Area, Step_Size, H,N )

%{
-------------------------   AUTHORSHIP  -------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

-------------------------   DESCRIPTION   -------------------------

This function creates a radio cell of diameter = Diameter and places
a maximum of 20 helpers (H) in well established positions in the cell.
N users are scattered randomly in the cell space with the condition that
any they cannot be closer than Radius_Protected_Area. The axis are
discretized to reduce the computation complexity.

------------------------- INPUT PARAMETERS -------------------------

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

------------------------- OUTPUT PARAMETERS -------------------------

-- Distances --
2D matrix of distances. The i-th row contains the distance of user i from
all the content providers in the cell. There are as many rows as the
number of users in the cell.

 
-- Fig_Handle -- 
A handle to a plot of the helpers, users and base station.

------------------------- EXAMPLE OF CALL -----------------------

A Radio Cell of 1 Km in diamter
10 meters of clear space around each user
Step size of 0.1 
20 Helpers in the network 
100 users

[ Distances ] = Place_Users( 1000, 10, 0.1, 20, 100)

%}


% ----------------------------   CODE     --------------------------

if H>20
    error('The system is not designed for more than 20 helpers.');
end
if (Step_Size > 1)
    error(['The script does not tollerate step sizes larger than 1.',...
           'If you wish to reduce the cell area, you can do so by changing the Diameter.']);
end

% Calculate the number of points of size Step_Size that can fit in Diameter
Points = Diameter/Step_Size;

% Initialize discrete axis 
x_axis = 1:1:Points;
y_axis = 1:1:Points;

% Define the centre of both axis
centre_x = round(Points/2);
centre_y = centre_x;

% Initialize the grid as a matrix of size Points x Points with value Inf.
% Any time a helper is placed in it, the value corresponding to its coordinates 
% is set to 1. For users is set to 0. This set of values (Inf, 1, 0) is completely 
% arbitrary and does not affect the behavior of the algorithm.

grid = Inf*ones(Points,Points);

% Initialize the matrix of distances
distances = zeros(N,H+1);   

% Calculate Rp, the radius of the protected area in units of Step_Size 
Rp = Radius_Protected_Area/Step_Size;


% Place Base station in the centre of the cell
bs_x = centre_x;
bs_y = centre_y;

% The grid where helpers are placed is obtained by dividing the space
% in 8 section vertically and 8 sections horizontally. We therefore define 
% the helpers unit 'hu' as 1/8 th of the total number of points and create a
% bidimentional array of 20 x,y coordinates. 

hu = floor(Points/8);

h_xy = zeros(20,2);
h_xy(1,:) = [3*hu,5*hu];
h_xy(2,:) = [5*hu,5*hu];
h_xy(3,:) = [5*hu,3*hu];
h_xy(4,:) = [3*hu,3*hu];
h_xy(5,:) = [2*hu,6*hu];
h_xy(6,:) = [4*hu,6*hu];
h_xy(7,:) = [6*hu,6*hu];
h_xy(8,:) = [6*hu,4*hu];
h_xy(9,:) = [6*hu,2*hu];
h_xy(10,:) = [4*hu,2*hu];
h_xy(11,:) = [2*hu,2*hu];
h_xy(12,:) = [2*hu,4*hu];
h_xy(13,:) = [1*hu,7*hu];
h_xy(14,:) = [4*hu,7*hu];
h_xy(15,:) = [7*hu,7*hu];
h_xy(16,:) = [7*hu,4*hu];
h_xy(17,:) = [7*hu,1*hu];
h_xy(18,:) = [4*hu,1*hu];
h_xy(19,:) = [1*hu,1*hu];
h_xy(20,:) = [1*hu,4*hu];


% Extract the coordinates of the first H helpers:
sel_helpers_xy = h_xy(1:H,:);

% Create an xy matrix for all providers (base station + helpers):
providers_xy = [sel_helpers_xy; bs_x, bs_y];

% Assign 1 to all those posiitons on the grid where there is provider:
for i = 1:1:H
   grid( sel_helpers_xy(i,1), sel_helpers_xy(i,2)) = 1;
end


%{
figure;
[hx,hy] = find(Grid == 1); 
plot(hx,hy,'o');
axis([0 Points 0 Points]);
%}


% Define uniform probability mass functions (pmf) for x and y
pmf_x = (1/Points)*ones(1,Points);
pmf_y = (1/Points)*ones(1,Points);

% Initialize seed to the rand function with the current clock time
rand('state', 100*sum(clock));

% Initialize matrix of coordinates for the N users.
placed_users_xy = zeros(N,2); 

% User index
i = 1;

% Index to track the number of iterations needed to place a user while
% skip == false (see below)
iteration = 1;

% Initialzie time 
time = 0;
avg_iter_per_placement = 1;
remaining_users = N;
while remaining_users > 0
    
    tic; % Timing
    
    % Randomly pick a user position on x and y according to the mass
    % probability functions Probability_x and Probability_y
    user_x = randsrc(1,1,[ x_axis; pmf_x ]);
    user_y = randsrc(1,1,[ y_axis; pmf_y ]);
    
    % Set a warning flag to true if the selected user is too close to the BS: 
    invalid_user_x = user_x > (centre_x - Rp ) & user_x < (centre_x + Rp );
    invalid_user_y = user_y > (centre_y - Rp ) & user_y < (centre_y + Rp );
    warning_bs_dist = invalid_user_x & invalid_user_y;
    
    % Set a warning flag to true if the choosen position for the user
    % overlaps with the positions assigned to helpers.
    warning_u_on_h = sum(sel_helpers_xy(:,1) == user_x  & sel_helpers_xy(:,2) == user_y) > 0;
    
    % Set a warning flag to true if the choosen position for the user
    % overlaps with an already placed user.
    warning_existing_u = sum(placed_users_xy(:,1) == user_x  & placed_users_xy(:,2) == user_y) > 0; 
    
    % Calculate the user-helpers distance:
    uh_delta_x = (sel_helpers_xy(:,1) - user_x).*Step_Size;
    uh_delta_y = (sel_helpers_xy(:,2) - user_y).*Step_Size;
    uh_distances = sqrt( abs( uh_delta_x ).^2 +  abs( uh_delta_y ).^2);
    
    % Set a warning flag to true if the placed user is too close to any of
    % the helpers:
    
    warning_h_dist = sum( uh_distances < Radius_Protected_Area );
    
    skip = warning_bs_dist | warning_u_on_h | warning_existing_u |warning_h_dist;
    if (skip)
        % The user is in the prohibited area. Keep selecting.
        i = i; % Skip.
    else  
        % Calculate the user-providers deltax and delta y
        up_delta_x = abs(providers_xy(:,1) - user_x);
        up_delta_y = abs(providers_xy(:,2) - user_y);
        
        % Caluclate Euclidean distances
        up_distances = sqrt(up_delta_x.^2 + up_delta_y.^2);
        
       % Make sure that Helpers can actually help this user. If the
       % distance between user and helpers is larger than that between user
       % and base-statio set the distance to Infinite (equivalent to
       % deactivate to helper becasue the SNR would be null). Never disable
       % the base station whose distance from the user is stored in the last
       % element of 'up_distances'.
       
        up_distances ( [ up_distances(1:end-1) >= up_distances(end); false] ) = Inf;
        
        % Solutions to the caching problem when users are connected only to one 
        % helper are triial. For this reason, the randomply placed user is discarded
        % if it is not connected to at least two helpers:
        
        user_gegree = sum(up_distances(1:end-1) < Inf);
        if ( user_gegree < 2 )
            i = i; % Skip.            
        else
            % save coordinates for this user!
            placed_users_xy(i,:) = [user_x, user_y];
            
            % Save the user-provider distances in the i-th row of the 'distances' matrix
            distances(i,:) = up_distances;
            
            % Set the grid posiiton corresponfing to where the user is to zero.       
            grid( user_x, user_y) = 0;
            avg_iter_per_placement = iteration/i;
            remaining_users = remaining_users - 1;
            i = i + 1;
        end
    end  
    time = time + toc;
    Avg_Time = time/iteration;
    Time_Left = avg_iter_per_placement * Avg_Time * remaining_users;
    disp(['USER PLACEMENT. Time to convergence: ' num2str(Time_Left) '.']);
    iteration = iteration + 1;
end

helpers_degree = sum(distances < Inf,1);


% PLOT THE NETWORK.

% Helpers are red circles
% Users are blue crosses
% Base station is a balck diamond

figure;
hold on;
title(['Random Users Placement with N=' num2str(N)  'users and H = ' num2str(H) ' helpers.']);
xlabel('Cell x_axis');
ylabel('Cell_y_axis');
axis([0 Points 0 Points]);
[hx,hy] = find(grid == 1); 
plot(hx,hy,'ro');
hold on;
[ux,uy] = find(grid == 0); 
plot(ux,uy,'x');
plot(bs_x, bs_y,'kd');

Fig_Handle = gcf;

end

