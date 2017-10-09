function [ Learning ] = Learners_Allocation( Learning_Setup, S, H, M, F )
% This function just allocates 

switch Learning_Setup
    
    case 'Single'    
        Initial_P = (1/F).*ones(1,F);   % Initial value probability mass function
        Initial_D = zeros(1,F);         % Initial estimates of the reward probability.
        Learner = struct('Search_Space',S, 'P',Initial_P, 'Ai',0, 'Z',zeros(1,F), 'W',zeros(1,F), 'D',Initial_D);
        for k = 1:1:M
            for j = 1:1:H
                Learning(k,j) = Learner;
            end
        end
    case 'Cooperative_Hard'
        
        % The H*M Files are HARD DIVIDED in M chunks of H files each.
        % Therefore, we have M search spaces, of smaller size. 
        Search_Spaces = zeros(M,H);   
        for k = 1:1:M
            st = 1 +(k-1)*H;
            ed = k*H;
            Search_Spaces(k,:) = S( st : ed );
        end
        Initial_P = (1/H).*ones(1,H);   % Initial value probability mass function
        Initial_D = zeros(1,H);         % Initial estimates of the reward probability.
        for k = 1:1:M
            for j = 1:1:H
                Learning(k,j) = struct('Search_Space', Search_Spaces(k,:), 'P',Initial_P, 'Ai',0, 'Z',zeros(1,H), 'W',zeros(1,H), 'D',Initial_D);
            end
        end
        
    case 'Cooperative_Shuffle'
        
        % The H*M Files are DIVIDED in M chunks of H files each RANDOMLY.
        % Therefore, we have M search spaces, of smaller size. 
        
        Initial_P = (1/H).*ones(1,H);   % Initial value probability mass function
        Initial_D = zeros(1,H);         % Initial estimates of the reward probability.


        for j = 1:1:H
            
            Actions_Left = S;
            for memory = 1:1:M
                for Nact = 1:1:H
                    Prob_Mass_Function = (1 / size(Actions_Left,2) )*ones(1, size(Actions_Left,2) );
                    Action = randsrc(1,1,[ Actions_Left; Prob_Mass_Function ]);
                    [~, cl] = find (Actions_Left == Action );
                    Actions_Left(cl) = [];
                    Search_Spaces(memory,Nact) = Action;
                end
            end
            
            for k = 1:1:M
                Learning(k,j) = struct('Search_Space', Search_Spaces(k,:), 'P',Initial_P, 'Ai',0, 'Z',zeros(1,H), 'W',zeros(1,H), 'D',Initial_D);
            end
        end

end


end

