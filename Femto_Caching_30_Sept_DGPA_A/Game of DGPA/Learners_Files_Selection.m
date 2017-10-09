function [ Available_Files, OUT_Learning] = Learners_Files_Selection( S, IN_Learning )

%% Learners Select Files in Parallel (same time)

OUT_Learning = IN_Learning;
NP = size(IN_Learning,2) + 1;
NL = size(IN_Learning,1);
NH = NP - 1;

Available_Files = zeros(NL,NH);

for j = 1:1: NH
    for k = 1:1:NL    
        Action = randsrc(1,1,[ S; IN_Learning(k,j).P ]);
        Available_Files(k,j) = Action;
        OUT_Learning(k,j).Ai = Action;
        OUT_Learning(k,j).Z(Action) = OUT_Learning(k,j).Z(Action) + 1;
    end
end

end

