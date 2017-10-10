function [ Available_Files, OUT_Learning] = Learners_Select_Files( IN_Learning )

%% Learners Select Files in Parallel (same time)

OUT_Learning = IN_Learning;
NP = size(IN_Learning,2) + 1;
NL = size(IN_Learning,1);
NH = NP - 1;

Available_Files = zeros(NL,NH);


for j = 1:1: NH
    for k = 1:1:NL
        Space = IN_Learning(k,j).Search_Space;
        PMF = IN_Learning(k,j).P ;
        Action = randsrc(1,1,[ Space; PMF ]);
        Available_Files(k,j) = Action;
        OUT_Learning(k,j).Ai = Action;
        [rw, cl] = find( Space == Action);
        OUT_Learning(k,j).Z(cl) = OUT_Learning(k,j).Z(cl) + 1;
    end
end
 

end

