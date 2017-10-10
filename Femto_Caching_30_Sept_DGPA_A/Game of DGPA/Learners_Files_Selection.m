function [ Available_Files, OUT_Learning] = Learners_Files_Selection( S, IN_Learning )

% Function that simulates the independent file selection by means of each learner. 

% IN_Learning is the matrix of Learners in input (before they select files)
% OUT_Learning is the matrix of Learners in output (after they select files)

OUT_Learning = IN_Learning;

% Number of providers
NP = size(IN_Learning,2) + 1;

% Number of learners per helper
NL = size(IN_Learning,1);

% Number of helpers
NH = NP - 1;

Available_Files = zeros(NL,NH);

for j = 1:1: NH
    for k = 1:1:NL   
        
        % Select file at random from library S.
        
        Action = randsrc(1,1,[ S; IN_Learning(k,j).P ]);
        
        Available_Files(k,j) = Action;
        
        % Save action in the learner 
        OUT_Learning(k,j).Ai = Action;
        OUT_Learning(k,j).Z(Action) = OUT_Learning(k,j).Z(Action) + 1;
        
    end
end

end

