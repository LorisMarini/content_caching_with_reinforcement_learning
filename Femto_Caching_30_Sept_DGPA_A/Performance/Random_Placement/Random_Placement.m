function [ Random_Allocation ] = Random_Placement( NL,NH )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

F = NL*NH;
S = 1:1:F;

Random_Allocation = zeros(NL,NH);

for j = 1:1: NH
    for k = 1:1:NL    
        Action = randsrc(1,1,[ S; (1/F)*ones(1,F) ]);
        Warning = sum( Random_Allocation(:,j) == Action) > 0;
        while Warning
            Action = randsrc(1,1,[ S; (1/F)*ones(1,F) ]);
            Warning = sum( Random_Allocation(:,j) == Action) > 0;
        end
        Random_Allocation(k,j) = Action;      
    end
end

end

