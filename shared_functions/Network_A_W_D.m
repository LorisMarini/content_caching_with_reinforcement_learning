function [ Net_A_W_D ] = Network_A_W_D( Network_Delays, S, Popularities, Allocation )
%{
-------------------------   AUTHORSHIP  -------------------------

Developer: Loris Marini
Affiliation: The University of Sydney
Contact: mrnlrs.tor@gmail.com
Notes:

-------------------------   DESCRIPTION   -------------------------

This function calculates the network average weighted delay as per eq(5) of
Marini, L., Li, J., & Li, Y. (n.d.). Distributed Caching based on Decentralized Learning Automata, 1–6.

------------------------- INPUT PARAMETERS -------------------------

-- Network_Delays --
Matrix of latencies [s] of size NxP with N = number of users, and P = number of 
content providers. Row i-th corresponds to user i, and each column to the download 
latency from corresponding provider (helper or base station).

-- S -- 


-- Popularities -- 

-- Allocation --


------------------------- OUTPUT PARAMETERS -------------------------

-- Net_A_W_D -- 
The average weighted delay of the entire network.

------------------------- EXAMPLE OF CALL -----------------------


% ----------------------------   CODE     --------------------------
%}

N = size(Network_Delays,1);
Users_W_Ds = zeros(1,N);

for n = 1:1:N
    User_Selections = User_NCA_Selection( n, S, Allocation, Network_Delays);
    Users_W_Ds(n) = User_Weighted_Delay( User_Selections, Popularities );
end

Net_A_W_D = sum(Users_W_Ds)/N;
end