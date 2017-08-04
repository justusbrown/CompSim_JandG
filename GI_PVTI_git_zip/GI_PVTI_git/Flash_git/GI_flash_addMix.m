function [success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,Zgas_vap, Zgas_liq, cubic_time]=GI_flash_addMix(components,Zi,Temp,pressure,thermo,options)

[components, T, p]=deal(components,Temp,pressure);
mixture=addmixture(components,T,p);
mixture.Zi=Zi;

%include numcells in initState and include it as input for flash
for i=1:numCells
    totalFluid{i}=addmixture(components, T, p);

[stability_flag_l,stability_flag_g,stability_flag, Zgas_vap, Zgas_liq] = stabilityTest(mixture,thermo);
success_flag = 1;
if stability_flag == 2
    success_flag = 0;
    vapor_frac=nan;
    vapor_y = nan; %
    liquid_x = nan;
    cubic_time=nan;
    Zgas_vap=nan;
    Zgas_liq=nan;
    
elseif stability_flag == 1 % which is single phase
    phase_flag = phase_Identify(mixture);
    if phase_flag ==1 %liquid single phase
        vapor_frac=0;
        vapor_y = zeros(size(mixture.Zi)); %
        liquid_x = mixture.Zi;
        Zgas_vap=nan;
        Zgas_liq=Zgas_liq;
        cubic_time=0;
    elseif phase_flag ==2 %vapor single phase
        vapor_frac=1;
        vapor_y = mixture.Zi; %
        liquid_x = zeros(size(mixture.Zi));
        Zgas_vap=Zgas_vap;
        Zgas_liq=nan;
        cubic_time=0;
    else %critical fluid
            success_flag = 0;  %needs to fix later
            vapor_frac=nan;
            vapor_y = nan; 
            liquid_x = nan;
            Zgas_vap=nan;
            Zgas_liq=nan;
            cubic_time=nan;
    end
else  % unstable
    [success_flag,vapor_y,liquid_x,vapor_frac,Zgas_vap, Zgas_liq, cubic_time]=vle2fash(mixture,thermo,options);
end



end