function [ Greedy_Allocation, Average_Weighted_Delays] = Greedy_Placement( Network_Delays, M, Popularities )
% Greedy Algorithm

H = size(Network_Delays,2) - 1;         % Number of Helpers.
N = size(Network_Delays,1);             % Number of Users.
S = 1:1:H*M;                            % Set of files we can off-load.
Available_Files = zeros(M,H);
N_Groups = nchoosek(S,M);
Min_Weighted_Delay = Inf*ones(1,N);
Remaining_Helpers = 1:1:H;
P_Selection = (1/H)*ones(1,H);
Total_Time = 0;


% Randomly select a helper.
Allocations_Left = size(N_Groups,1)*H;
Counter = 1;
Average_Weighted_Delays = zeros(1,H);
Weighted_Delays = zeros(1,N);
    
while size(Remaining_Helpers,2) ~= 0
    
    j = randsrc(1,1,[ Remaining_Helpers; P_Selection ]);
    Remaining_Helpers(Remaining_Helpers == j) = [];    
    P_Selection = (1/size(Remaining_Helpers,2) )*ones(1, size(Remaining_Helpers,2) );
    for i = 1:1: size(N_Groups,1)
        tic;
        Available_Files(:,j) = N_Groups(i,:); 
        for n = 1:1:N
            %tic;
            User_Selections = User_NCA_Selection( n, S, Available_Files, Network_Delays);
            %toc
            %tic;
            Weighted_Delays(n) = User_Weighted_Delay( User_Selections, Popularities );
            %toc
        end    
        Current_Helper_Degree = sum( Network_Delays(:,j) < Inf );
        Average_Weighted_Delays(j) = sum( Weighted_Delays( Network_Delays(:,j) < Inf ) ) / Current_Helper_Degree;
        if ( Average_Weighted_Delays(j) < Min_Weighted_Delay )
            Min_Weighted_Delay = Average_Weighted_Delays(j);
            Current_Best_Allocation = N_Groups(i,:);
        end
        Allocations_Left = Allocations_Left - 1;
        Total_Time = Total_Time + toc;
        Avg_Partial_Time = Total_Time/Counter;
        Time_Left = Avg_Partial_Time * Allocations_Left;
        Counter = Counter +1;
        disp(['GREEDY Placement. Time left: ' num2str(floor(Time_Left/60)), ' minutes and ' num2str(rem(Time_Left,60)) ' seconds.']);
    end
    Available_Files(:,j) = Current_Best_Allocation;    
end

Greedy_Allocation = Available_Files;

end

