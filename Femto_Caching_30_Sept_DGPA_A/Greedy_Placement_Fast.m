function [ Available_Files, Helpers_Avg_Weighted_Delay ] = Greedy_Placement_Fast( Network_Delays, M, Popularities )

% GREEDY ALGORITHM OPTIMISED FOR SPEED.

H = size(Network_Delays,2) - 1;         % Number of Helpers.
N = size(Network_Delays,1);             % Number of Users.
F = M*H;
S = 1:1:F;                              % Set of files we can off-load.
Available_Files = zeros(M,H);
N_Groups = nchoosek(S,M);
Remaining_Helpers = 1:1:H;
P_Selection = (1/H)*ones(1,H);
Allocations_Left = size(N_Groups,1)*H;
Weights = repmat(Popularities, N, 1);
Helpers_Avg_Weighted_Delay = Inf*ones(1,H);

Counter = 1;
Total_Time = 0;

while size(Remaining_Helpers,2) ~= 0
    
    j = randsrc(1,1,[ Remaining_Helpers; P_Selection ]);
    Remaining_Helpers(Remaining_Helpers == j) = [];
    P_Selection = (1/size(Remaining_Helpers,2) )*ones(1, size(Remaining_Helpers,2) );
    Min_Weighted_Delay = Inf;
    
    for i = 1:1: size(N_Groups,1)
        tic;
        Available_Files(:,j) = N_Groups(i,:);
        Users_NCA_Delays = zeros(N,F);
        for f = 1:1:F
            Altern_From_Helpers = (Available_Files == f);              
            Who_Has_Not_This_File = sum(Altern_From_Helpers,1) == 0;
            Delays = Network_Delays;
            Delays(:, [Who_Has_Not_This_File, false] ) = Inf;
            Users_NCA_Delays(:,f) = min(Delays,[],2);
        end
        Weighted_Delays = sum(Users_NCA_Delays .* Weights,2);
        Current_Helper_Degree = sum( Network_Delays(:,j) < Inf );
        Helpers_Avg_Weighted_Delay(j) = sum( Weighted_Delays( Network_Delays(:,j) < Inf ) ) / Current_Helper_Degree;
        
        if ( Helpers_Avg_Weighted_Delay(j) < Min_Weighted_Delay )
            Min_Weighted_Delay = Helpers_Avg_Weighted_Delay(j);
            Current_Best_Allocation = N_Groups(i,:);
        end
        Allocations_Left = Allocations_Left - 1;
        Total_Time = Total_Time + toc;
        Avg_Partial_Time = Total_Time/Counter;
        Time_Left = Avg_Partial_Time * Allocations_Left;
        Counter = Counter +1;
        %disp(['GREEDY Placement. Time left: ' num2str(floor(Time_Left/60)), ' minutes and ' num2str(rem(Time_Left,60)) ' seconds.']);
    end
    Available_Files(:,j) = Current_Best_Allocation;
end
end



