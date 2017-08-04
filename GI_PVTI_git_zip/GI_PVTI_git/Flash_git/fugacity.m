function [liq_fug, vap_fug,Zgas_liq, Zgas_vap, fug_flag] = fugacity(mixture, thermo, liquid_x, vapor_y)

fug_flag=1;
thermo.fugacity_switch=1; % switch on
eosf = thermo.EOS;

%------------------------------------
%liquid fugacity
thermo.phase=1;
mixture.Zi = liquid_x;
[liq_fug_coef,success_flag_l, Zgas_liq]=eosf(mixture, thermo);

%---------------------------------------
%vapor fugacity
thermo.phase=2;
mixture.Zi = vapor_y;
[vap_fug_coef,success_flag_g,Zgas_vap]=eosf(mixture, thermo);
%--------------------------------------
if success_flag_l ==0 || success_flag_g ==0
    fug_flag = 0;
    liq_fug = nan;
    vap_fug = nan;
else
    p = mixture.pressure; %[Pa]
    liq_fug = liq_fug_coef.*liquid_x*p; %[Pa]
    vap_fug = vap_fug_coef.* vapor_y*p; %[Pa]
    Zgas_liq=Zgas_liq;
    Zgas_vap=Zgas_vap;
end

end

