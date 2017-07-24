%%
%
% perlOneLineDescription(Returns quadratic relative permeability)
% 
%%


function [krL, krG] = quadraticRelPerm_JandG(So)
    krL = So.*So;
    krG = (1 - So).^2;
 end