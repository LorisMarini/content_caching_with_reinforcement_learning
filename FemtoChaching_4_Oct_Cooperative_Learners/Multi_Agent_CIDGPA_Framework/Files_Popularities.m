function [ Popularities ] = Files_Popularities( S, Gamma_ZipF )

% Gamma_ZipF: Determination of Popularities: ZIPF Distribution Law.

 Popularities = 1./(S.^Gamma_ZipF)...  % We assume that S corresponds to the action's ranking. Action 1 is then
     ./(sum((1./S).^Gamma_ZipF));      % the most popular file. Action F is the least popular file.
end

