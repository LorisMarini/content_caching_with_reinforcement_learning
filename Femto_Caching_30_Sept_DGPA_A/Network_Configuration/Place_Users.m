function [ Distances, Fig_Handle ] = Place_Users( Diameter,Radius_Protected_Area, Step_Size, H,N )
% PLACE USERS RANDOMLY

% H = Number of Helpers;
% N = Number of users;
% Diameter = Diameter (edge of the square) in meters.
% Step_Size =     In meters always lesser or equal to 1.

% EXAMPLE OF CALL:
% [ Distances ] = Users_Placement( 1000, 10, 0.1, 20, 100)

if H>20
    error('The system is not designed for more than 20 helpers.');
end
if (Step_Size > 1)
    error('The script does not tollerate step sizes larger than 1. If you wish to reduce the cell area, you can do so by changing the Diameter.');
end
Distances = zeros(N,H+1);             
Length_Radius_Protected_Area = Radius_Protected_Area/Step_Size;
R = Length_Radius_Protected_Area;
Points = Diameter/Step_Size;
x_axis = 1:1:Points;
y_axis = 1:1:Points;
Probability_x = (1/Points)*ones(1,Points);
Probability_y = (1/Points)*ones(1,Points);
Median_x = round(Points/2);
Median_y = Median_x;

Grid = Inf*ones(Points,Points);

BS_x = Median_x;
BS_y = Median_y;

u = floor(Points/8);
Helpers_Grid_Coordinates = zeros(20,2);
Helpers_Grid_Coordinates(1,:) = [3*u,5*u];
Helpers_Grid_Coordinates(2,:) = [5*u,5*u];
Helpers_Grid_Coordinates(3,:) = [5*u,3*u];
Helpers_Grid_Coordinates(4,:) = [3*u,3*u];
Helpers_Grid_Coordinates(5,:) = [2*u,6*u];
Helpers_Grid_Coordinates(6,:) = [4*u,6*u];
Helpers_Grid_Coordinates(7,:) = [6*u,6*u];
Helpers_Grid_Coordinates(8,:) = [6*u,4*u];
Helpers_Grid_Coordinates(9,:) = [6*u,2*u];
Helpers_Grid_Coordinates(10,:) = [4*u,2*u];
Helpers_Grid_Coordinates(11,:) = [2*u,2*u];
Helpers_Grid_Coordinates(12,:) = [2*u,4*u];
Helpers_Grid_Coordinates(13,:) = [1*u,7*u];
Helpers_Grid_Coordinates(14,:) = [4*u,7*u];
Helpers_Grid_Coordinates(15,:) = [7*u,7*u];
Helpers_Grid_Coordinates(16,:) = [7*u,4*u];
Helpers_Grid_Coordinates(17,:) = [7*u,1*u];
Helpers_Grid_Coordinates(18,:) = [4*u,1*u];
Helpers_Grid_Coordinates(19,:) = [1*u,1*u];
Helpers_Grid_Coordinates(20,:) = [1*u,4*u];

Helpers = Helpers_Grid_Coordinates(1:H,:); % Extract only the Helpers we are interested in.
Providers = [Helpers; BS_x, BS_y];

for i = 1:1:H
   Grid( Helpers(i,1), Helpers(i,2)) = 1;
end


%{
figure;
[hx,hy] = find(Grid == 1); 
plot(hx,hy,'o');
axis([0 Points 0 Points]);
%}

rand('state', 100*sum(clock));
Placed_Users = zeros(N,2); % Coordinates of the placed users.
Users_Left = N;
i = 1;
time = 0;
iteration = 0;
Avg_Iterations_per_placement = Inf;

while Users_Left > 0
    tic;
    user_x = randsrc(1,1,[ x_axis; Probability_x ]);
    user_y = randsrc(1,1,[ y_axis; Probability_y ]);
    % We make sure that the user selected is not in the proximity of the BS: 
    Warning_BS = user_x > (Median_x - R ) & user_x < (Median_x + R ) & user_y > (Median_y - R ) & user_y < (Median_y + R );
    Warning_H = sum(Helpers(:,1) == user_x  & Helpers(:,2) == user_y) > 0;
    Warning_U = sum(Placed_Users(:,1) == user_x  & Placed_Users(:,2) == user_y) > 0; 
    
    User_Distances_WithRespectToThe_Helpers = sqrt( abs( (Helpers(:,1) - user_x).*Step_Size ).^2 +  abs( (Helpers(:,2) - user_y).*Step_Size ).^2);
    Warning_H_Dist = sum( User_Distances_WithRespectToThe_Helpers < Radius_Protected_Area );
    
    Skip = Warning_BS | Warning_H | Warning_U |Warning_H_Dist;
    if (Skip)
        % The user is in the prohibited area. Keep selecting.
        i = i; % Skip.
    else  
        Delta_x = abs(Providers(:,1) - user_x);
        Delat_y = abs(Providers(:,2) - user_y);
        Dist = sqrt(Delta_x.^2 + Delat_y.^2);
       % Let's make sure that Helpers can actually help.
        Who_To_Prune = [ Dist(1:end-1) >= Dist(end); false];
        Dist (Who_To_Prune) = Inf;
        
        % If the user is not connected to at least two helpers, it should be skipped.
        User_Degree = sum(Dist(1:end-1) < Inf);
        if ( User_Degree < 2 )
            i = i; % Skip.            
        else
            Placed_Users(i,:) = [user_x, user_y];
            Distances(i,:) = Dist;
            Grid( user_x, user_y) = 0;
            Avg_Iterations_per_placement = (iteration)/i;
            Users_Left = Users_Left - 1;
            i = i +1;
        end
    end  
    time = time + toc;
    Avg_Time = time/iteration;
    Time_Left = Avg_Iterations_per_placement*Avg_Time*Users_Left;
    disp(['USER PLACEMENT. Time to convergence: ' num2str(Time_Left) '.']);
    iteration = iteration + 1;
end

Helpers_Degree = sum(Distances < Inf,1);
%{
if ( sum(Helpers_Degree < 2) > 1 )
    error('We envisioned a minimum number of 2 users connected to each helper.');
end
%}

% PLOT THE NETWORK.
figure;
hold on;
title(['Random Users Placement with N=' num2str(N)  'users and H = ' num2str(H) ' helpers.']);
xlabel('Cell x_axis');
ylabel('Cell_y_axis');
axis([0 Points 0 Points]);
[hx,hy] = find(Grid == 1); 
plot(hx,hy,'ro');
hold on;
[ux,uy] = find(Grid == 0); 
plot(ux,uy,'x');
plot(BS_x, BS_y,'kd');

Fig_Handle = gcf;

end

