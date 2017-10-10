function [ Convergence ] = Check_Game_Convergence( Learning, P_Threshold )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Learners, Helpers] = size(Learning);
C_Flag = zeros(Learners, Helpers);

for j=1:1:Helpers
    for k=1:1:Learners
        
        if sum(Learning(k,j).P > P_Threshold) >= 1
            C_Flag(k,j) = 1;
        else
            C_Flag(k,j) = 0;
        end
    end
end

if (sum(sum(C_Flag == 1)) == Helpers*Learners)
    
    Convergence = 1;
else
    Convergence = 0;
end

end

