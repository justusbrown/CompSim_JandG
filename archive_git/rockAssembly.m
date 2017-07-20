mrstModule add ad-fi 
mrstModule add ad-core


G = cartGrid([2, 2, 2]);
G = computeGeometry(G);
% Disable gravity
gravity off
% Set up uniform permeability and constant porosity
rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = [0.5, 0.1,0.3,0.4,0.8,0.7,0.2,0.3]';
rock.pv=poreVolume(G,rock);
T=computeTrans(G,rock);

%%
%SETUP THE CONTROLS
%OPTIONS FOR PVT
options.convergence_eps = 1e-12;   %convergence tolerance for fugacity
options.trivial_eps = 1e-3;     %trivial shift for bisection algorithm
options.RRiteration = 200;   %maximum number of Rachford Rice iteration using Newton's method
options.max_outer_loop = 1000;   %max number of fugacity updates

influx_p=13e6;
outflux_p=8e6;

influx_rate = 1000/day  % in m^3/s


Temp = 30 + 273.15;  % Temperature (in Kelvin)
thermo=addThermo();
thermo.EOS=@PREOS;


influx_C  = [2.12,0.06,0.56,4.43];
influx_C  =mynormalize(influx_C);
outflux_C = [2.12,0.06,0.56,4.43];
outflux_C=mynormalize(outflux_C);
C_initial = [2.12,0.06,0.56,4.43]; 
C_initial=mynormalize(C_initial);
[component, comp_flag] = addComponents({'CO2','CH4','C10H22','H2O'})

influxFluid=addMixture(component,Temp,influx_p);
influxFluid.mole_fraction=influx_C;
outfluxFluid=addMixture(component,Temp,outflux_p);
outfluxFluid.mole_fraction=outflux_C;
initialFluid=addMixture(component,Temp,outflux_p);
initialFluid.mole_fraction=C_initial;
clf
plotCellData(G,G.cells.indexMap), view(30,50), xlabel('x'), ylabel('y'),zlabel('z'), colorbar

%%
%CONFUSED ON THIS BLOCK. COULD JUST SAY bc.dirichlet.faces=2
bound_cell_out=1;
bound_faces = G.cells.faces(G.cells.facePos(bound_cell_out):G.cells.facePos(bound_cell_out + 1) - ...
                              1 , :);
bc.dirichlet.faces =  bound_faces(bound_faces(:, 2) == 2 , 1);
%%
%%
%SET UP DIRICHLET CONDITION
bc.dirichlet.pressure = outflux_p;
bc.dirichlet.fluid=outfluxFluid
% flash for surface properties
[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,cubic_time]=GI_flash(bc.dirichlet.fluid,thermo,options)
bc.dirichlet.Xig=vapor_y(1:3);
bc.dirichlet.Xio=liquid_x(1:3);
bc.dirichlet.So=1-vapor_frac;
bc.dirichlet.Zi=bc.dirichlet.Xig*vapor_frac+bc.dirichlet.Xio*bc.dirichlet.So;%NEED TO KNOW WHAT vapor_frac IS. (does it include water?)
%bc.dirichlet.Eio=bc.dirichlet.Xio./rock.pv; %THIS IS SUPPOSED TO BE DEFINITION OF Eio, MOLAR DENSITY, BUT I DO NOT THINK IT WORKS
%bc.dirichlet.Eig=bc.dirichlet.Xig./rock.pv; %ALSO rock.pv IS A VECTOR CONTAINING PORE VOL OF EACH CELL, SO THIS NEEDS CHANGING, BUT NOT YET FOR CONCEP
bc.dirichlet.Eo=sum(bc.dirichlet.Xio)/rock.pv(1);
bc.dirichlet.Eg=sum(bc.dirichlet.Xig)/rock.pv(1); %HERE I AM TREATING MOLES AND MOLE FRACTIONS AS EQUIVALENT BUT I KNOW THATS NOT OKAY
bc.dirichlet.F=bc.dirichlet.Eo*bc.dirichlet.So+bc.dirichlet.Eg*vapor_frac; %STILL NEED TO CLARIFY So, Sg,Sw, vs L and V
bc.dirichlet.Sw=vapor_y(4)/sum(vapor_y)+liquid_x(4)/sum(liquid_x);

%bc.dirichlet for outflux IS NOW DEFINED FOR Zi, F, Sw, and P ... (OUR PRIMARY VARIABLES)
%NOW DEFINE DIRICHLET CONDITIONS FOR INFLUX
bound_cell_in=8;
bc.influx_cells=bound_cell_in;
[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,cubic_time]=GI_flash(influxFluid,thermo,options) %DONT SEEM TO USE ANYTHING FROM THIS FLASH TEST...WHY?
for ic=1:3 %DOING WATER SEPERATELY
    bc.C_influx(ic)=influx_rate*(vapor_y(ic)*vapor_frac+liquid_x(ic)*(1-vapor_frac)); %IM NOT SURE IF THIS IS CORRECT EITHER
end
bc.water_influx=influx_rate*(liquid_x(4)*(1-vapor_frac)+vapor_y(4)*vapor_frac);

%THIS IS THE END OF SETTING UP THE CONTROLS!!!

%%NOW NEED TO SETUP SYSTEM AND NONLINEAR SOLVER


%%
%NOW INITIALIZING THE STATE
[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,cubic_time]=GI_flash(initialFluid,thermo,options);

numCells=G.cells.num;

state.C=C_initial(1:4).*ones(numCells,1);
state.Xig=vapor_y(1:3).*ones(numCells,1);;
state.Xio=liquid_x(1:3).*ones(numCells,1);;
state.So=(1-vapor_frac).*ones(numCells,1);;
state.Zi=(state.Xig.*vapor_frac+state.Xio.*state.So).*ones(numCells,1);
state.Eo=sum(state.Xio)./rock.pv;
state.Eg=sum(state.Xig)./rock.pv;%HERE I AM TREATING MOLES AND MOLE FRACTIONS AS EQUIVALENT BUT I KNOW THATS NOT OKAY
state.Ew=(liquid_x(4)+vapor_y(4)).*ones(numCells,1)./rock.pv;
state.F=state.Eg.*vapor_frac+state.Eo.*state.So;
state.Sw=vapor_y(4)/sum(vapor_y)+liquid_x(4)/sum(liquid_x);
state.pressure=outflux_p*ones(numCells,1);

state0=state;

%THIS IS THE END OF INITIALIZING THE STATE
%%LETS START SETTING UP AND SOLVING THE SYSTEM
dt = 200*day;              % Time step
total_time = 100000*day;   % Total time

steps      = dt*ones(floor(total_time/dt), 1); 
t          = cumsum(steps); 


%% 
% HERE ARE THE TIME STEPS STARTIN FOR THE ITERATIONS

for tstep = 1 : numel(steps)

   dt = steps(tstep); 
   
   %THE MAIN VARIABLESARE P,F, Zi, Sw So MAKE ADI VARIABLES
   
   p=state.pressure;
   F=state.F;
   Zi=state.Zi;
   Sw=state.Sw;
   Xig=state.Xig;
   Xio=state.Xio;
   Ew=state.Ew;
   So=state.So;

   [p,F,Zi,Sw]=initVariablesADI(p,F,Zi,Sw);
   
   p0=state0.pressure;
   F0=state0.F;
   Zi0=state0.Zi
   Sw0=state0.Sw
   Xig0=state0.Xig;
   Xio0=state0.Xio;
   Ew0=state0.Ew;
   So0=state0.So;
   
[krL,krG]=quadraticRelPerm(So);
bd=bc.dirichlet;
[bc_krL, bc_krG] = quadraticRelperm(bd.So);





