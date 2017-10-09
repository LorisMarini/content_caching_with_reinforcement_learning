function [ Net_A_W_D ] = Network_A_W_D( Network_Delays, S, Popularities, Allocation )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

N = size(Network_Delays,1);
Users_W_Ds = zeros(1,N);

for n = 1:1:N
    User_Selections = User_NCA_Selection( n, S, Allocation, Network_Delays);
    Users_W_Ds(n) = User_Weighted_Delay( User_Selections, Popularities );
end

Net_A_W_D = sum(Users_W_Ds)/N;
end