function [ Learning ] = Allocate_Learners( H,M,F )

% This function just allocates MxH learners in a matrix. Each element of
% the matrix is a struct containing variables such as te initial
% probability vector P, The currently chosen action Ai, counters for the
% numebr of times actions are selected Z and rewarded W, and the vector
% containing the estimate sof the reward probabilities D.

Initial_P = (1/F).*ones(1,F);   % Initial value probability mass function
Initial_D = zeros(1,F);         % Initial estimates of the reward probability.

Learner = struct('P',Initial_P, 'Ai',0, 'Z',zeros(1,F), 'W',zeros(1,F), 'D',Initial_D);
 for k = 1:1:M
     for j = 1:1:H
         Learning(k,j) = Learner;
     end
 end

end

