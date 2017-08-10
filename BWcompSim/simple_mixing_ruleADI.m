function [a, b] = simple_mixing_ruleADI(mixture, thermo, ai, bi)
% 

bip = mixture.bip;
mixing_rule_num = thermo.mixingrule;
temperature = mixture.Temp;
%{
if iscell(mixture.Zi)
x=cell2mat(mixture.Zi.val);
else
%}
x = mixture.Zi;


bipeos = [bip.EOScons]+[bip.EOStdep]*temperature;
N = length(x.val);
if (mixing_rule_num == 1)  %simple van der Waals mixing
    b = x.val*bi';
    a=0;
       for i=1:N
         for j=1:N
           a=a+x(i)*x(j)*sqrt(ai(i)*ai(j))*(1-bipeos(i,j));
         end
       end
end

end
    
