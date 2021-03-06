function phase_flag = phase_Identify(mixture)



critical_temp = [mixture.components.Tc]; %[K]
z = mixture.Zi;
T = mixture.Temp;
Tc = z*critical_temp';
if T<Tc
    phase_flag = 1; %liquid
else
    phase_flag = 2; %vapor
end

end
