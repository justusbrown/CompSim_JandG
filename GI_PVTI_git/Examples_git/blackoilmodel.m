% test case
% liquid-vapor phase equilibia calculation
% Define the components and load pure physsical properties
z=[2.12,0.06,41.56,20.43,10.94,2.81,2.42,1.71,1.45,3.11,30.39];
z0=z
z=mynormalize(z);
[component, comp_flag] = addComponents({'CO2','N2','CH4','C2H6','C3H8','C4H10','C4H10','C5H12','C5H12','C6H14','C8H18'});
% Define the thermodynamic models
T= 360; % [K]  temperature
p= 20e6; % [Pa] pressure
thermo = addThermo();
thermo.EOS = @PREOS;
mixture = addMixture(component, T, p);
mixture.mole_fraction = z;
% Define flash options
options.convergence_eps = 1e-12;   %convergence tolerance for fugacity
options.trivial_eps = 1e-3;     %trivial shift for bisection algorithm
options.RRiteration = 200;   %maximum number of Rachford Rice iteration using Newton's method
options.max_outer_loop = 1000;   %max number of fugacity updates
% flash for surface properties
[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,Zgas_vap, Zgas_liq, cubic_time]=GI_flash(mixture,thermo,options);

