function [ Weighted_Delay ] = User_Weighted_Delay( User_Selections, Popularities )

% Calculates the user weighted delay.

Weighted_Delay = 0;   % Average delay for user 'n'.
N_Files_Selected = 0; % Error Check.

NP = size(User_Selections,2); % Number of Providers (H+1)
NL = size(User_Selections,1); % Number of Learners
F = size(User_Selections,3);  % Number of Files

for f = 1:1:F
    FLG = 1;
    for j = 1:1: NP
        for k = 1:1: NL
            if ( j < NP)
                % When user 'n' has selected file f from a lerner in one helper
                if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                    Weighted_Delay = Weighted_Delay + Popularities(f).* User_Selections(k,j,f);
                    N_Files_Selected = N_Files_Selected +1;
                end
            elseif ( j == NP && FLG)
                % When user 'n' has selected file f from the Base Station
                if (User_Selections(k,j,f) ~= Inf && User_Selections(k,j,f) > 0)
                    Weighted_Delay = Weighted_Delay + Popularities(f).* User_Selections(k,j,f);
                    N_Files_Selected = N_Files_Selected +1;
                    FLG = 0;
                end
            end
        end % For all Lerners (k)
    end % For all Sources (j)
end % For all F files (f)

if (N_Files_Selected ~= F)
    error(['User ' num2str(n) ' did not select all the files in S when evaluating his weighted delay.']);
end
end

