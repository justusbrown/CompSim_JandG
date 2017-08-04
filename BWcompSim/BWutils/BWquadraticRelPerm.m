%%
%
% perlOneLineDescription(Returns quadratic relative permeability)
% 
%%


function [krL, krG] = BWquadraticRelPerm(So)
    krL = So.*So;
    krG = (1 - So).^2;
 end