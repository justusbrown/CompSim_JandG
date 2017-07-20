%I COPY AND PASTED BUT HAVENT DONE ANYTHING YET
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
rock.T=computeTrans(G,rock);
rock.G=G

%%
%SETUP THE CONTROLS
%OPTIONS FOR PVT
options.convergence_eps = 1e-12;   %convergence tolerance for fugacity
options.trivial_eps = 1e-3;     %trivial shift for bisection algorithm
options.RRiteration = 200;   %maximum number of Rachford Rice iteration using Newton's method
options.max_outer_loop = 1000;   %max number of fugacity updates

influx_p=10e6;
outflux_p=8e6;

influx_rate = 1000/day  % in m^3/s


Temp = 30 + 273.15;  % Temperaclearture (in Kelvin)
thermo=addThermo();
thermo.EOS=@PREOS;


influx_C  = [2,4,3,1];
influx_C  =mynormalize(influx_C);
outflux_C = [2,4,3,1];
outflux_C=mynormalize(outflux_C);
C_initial = [2,4,3,1]; 
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
%THIS IS WHERE setupControls WAS gr-07/19
bc=setupControls_JandG(rock,outfluxFluid,influxFluid,influx_rate,thermo,options);

%%
%SETUP SYSTEM
   cf = rock.G.cells.faces(:,1);
   nf = rock.G.faces.num;
   nc = rock.G.cells.num;
   %NOTE using rock.pv and rock.poro, don't want an s object like bravodome
   %COMPUTE FULLTRANSMITIBILITY
   rock.Tfull=1 ./ accumarray(cf,1./rock.T,[nf,1]);
   rock.T=rock.Tfull;
   
   %SETUPDISCRETE DIVERGENCE OPERATOR
   [C1,C2,C_div,div]=divOp_JandG(rock.G, cf, nf, nc);
   
   %GAGE 7/19 STARTED SETTING UP SYSTEM. FINSHED DIVERGENCE OPERATOR. 
   %CODE WORKS UP TO THIS POINT. NEED GRADIENT OPERATOR STILL.
   %DECIDED TO KEEP SETUP SYSTEM IN THE CORE FILE. THIS WILL NOT BE THE
   %CASEWITH ASSEMBLING THE EQUATIONS.
   %NEXT STEP IS SETTING UP GRADIENT OPERATOR
   
      % Set up the discrete gradient operator, |grad|. We compute the differences of cell values
   % across each face. It is a linear mapping from cells' to faces' values.
   N = double(rock.G.faces.neighbors);
   index = 1:nf;
   interior = prod(N, 2)~=0;

   C1_interior = sparse(index(interior), N(interior, 1), ones(nnz(interior), 1), nf, nc);
   C2_interior = sparse(index(interior), N(interior, 2), ones(nnz(interior), 1), nf, nc);

   
   %%
   % Compute the boundary contribution to the gradient operator. They corresponds to the
   % external faces where Dirichlet conditions are given. We are careful to use the
   % correct signs.
   % 

   is_dirichlet_faces1 = N(bc.dirichlet.faces, 1) ~= 0;
   is_dirichlet_faces2 = ~is_dirichlet_faces1;

   dirichlet_faces1 = bc.dirichlet.faces(is_dirichlet_faces1);
   dirichlet_faces2 = bc.dirichlet.faces(is_dirichlet_faces2);

   C1_exterior = sparse(index(dirichlet_faces1), ...
                        N(dirichlet_faces1, 1), ...
                        ones(numel(dirichlet_faces1), 1), nf, nc);
   C2_exterior = sparse(index(dirichlet_faces2), ...
                        N(dirichlet_faces2, 2), ...
                        ones(numel(dirichlet_faces2), 1), nf, nc);
                    
   %%
   % The gradient operator is the sum of internal and boundary contributions.
   %
   
   C = C1_interior + C1_exterior - (C2_interior + C2_exterior);

   pressure_bc = sparse(nf, 1);
   pressure_bc(dirichlet_faces1) = - bc.dirichlet.pressure(is_dirichlet_faces1);
   pressure_bc(dirichlet_faces2) = + bc.dirichlet.pressure(is_dirichlet_faces2);

   p_grad = @(p)(C*p + pressure_bc);  

   grad = @(val, bc_val)(grad_JandG(val, bc_val, nf, C, ...
                                is_dirichlet_faces1, dirichlet_faces1, ... 
                                is_dirichlet_faces2, dirichlet_faces2));

%%
   % Set up the gravity term.
   %
   
   z = rock.G.cells.centroids(:, 3);
   fz = rock.G.faces.centroids(:, 3);
   dz = grad(z, fz);
   
      %% 
   % Define function to compute concentrations on faces using cells' values and upwind
   % directions. Here, the _nature_ of the boundary faces (Dirichlet or note) are not
   % allowed to change, although the _values_ on these faces may change.
   %
   
   faceConcentrations = @(flag, conc_c, bc_conc) ...
       faceConcentrations_JandG(flag, conc_c, bc_conc, N, interior, dirichlet_faces2, ...
                          dirichlet_faces1,  bc, nf, nc);
  %THIS IS THE END OF SETUP SYSTEM. I BELIEVE IT IS COMPLETE, and KNOW 
  %THAT EVERYTHING STILL FUNCTIONS.
  
  %%NONLINEAR SOLVER PARAMETERS START HERE:
   %SETUP NONLINEAR SOLVER WILL START HERE
   maxIteration=50;
   nonlinear=setNonlinearSolverParameters_JandG(maxIteration);
   %Done setting up nonlinear solver gr 07/20
   
   %%NOW WE WILL INITIALIZE THE STATE USING initState_JandG
   %RIGHT AFTER INIT STATE IS BEING PUT IN NOW gr 07/20
   [state0]=initState_JandG(rock,initialFluid,options,thermo);




%{

%%
%NOW INITIALIZING THE STATE
[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,Zgas_vap, Zgas_liq, cubic_time]=GI_flash(initialFluid,thermo,options);

numCells=G.cells.num;

%FUTURE THOUGHTS:
%1) THIS SATURATIONS WILL NEED TO BE BACK CALCULATED FROM THE PRESENT F
%2) EPSILON IS EQUAL TO A RATIO OF MOLS/VOL OBTAINED FROM PREOS
%3) V AND L WILL ALSO BE FROM PREOS
%4) SW CAN THEN BE COMPUTED FROM OTHER SATURATIONS

state.C=C_initial(1:4).*ones(numCells,1);
state.Xig=vapor_y(1:3).*ones(numCells,1);
state.Xio=liquid_x(1:3).*ones(numCells,1);
state.So=(1-vapor_frac).*ones(numCells,1);
state.Zi=(state.Xig.*vapor_frac+state.Xio.*state.So).*ones(numCells,1);
state.Eo=sum(state.Xio)./rock.pv;
state.Eg=sum(state.Xig)./rock.pv;%HERE I AM TREATING MOLES AND MOLE FRACTIONS AS EQUIVALENT BUT I KNOW THATS NOT OKAY
state.Ew=(liquid_x(4)+vapor_y(4)).*ones(numCells,1)./rock.pv;
state.F=state.Eg.*vapor_frac+state.Eo.*state.So;
state.Sw=vapor_y(4)/sum(vapor_y)+liquid_x(4)/sum(liquid_x);
state.pressure=outflux_p*ones(numCells,1);
state.cwL=liquid_x(4).*ones(numCells,1);
state.cwV=vapor_y(4).*ones(numCells,1);
state.Cw=(vapor_y(4)*vapor_frac+liquid_x(4)*(1-vapor_frac)).*ones(numCells,1); %IM NOT SURE THIS IS ASSEMBLED CORRECTLY, BUT IM TRYING NOT TO FIXATE ON INDIVIDUAL LINES
state.Sg=vapor_frac;%NEED TO KNOW IF VAPOR_FRAC INCLUDES WATER. I THINK IT DOES AND SO SOME THINGS NEED FIXING

state.fluid=initialFluid
%CONSTANT VISCOSITIES & COMPRESSIBILITY (which I don't know if I need atm)
state.fluid.muL=1e-3;
state.fluid.muG=1e-5;
state.fluid.cl    = 4.4e-5/atm  % Compressibility
state.fluid.p_ref = 1*atm       % Reference pressure


state0=state;

%THIS IS THE END OF INITIALIZING THE STATE
%%LETS START SETTING UP AND SOLVING THE SYSTEM
dt = 200*day;              % Time step
total_time = 100000*day;   % Total time

steps      = dt*ones(floor(total_time/dt), 1); 
t          = cumsum(steps); 

nComp_C=3; %# OF NON WATER COMPONENTS
%% 
% HERE ARE THE TIME STEPS STARTINg FOR THE ITERATIONS

for tstep = 1 : numel(steps)

   dt = steps(tstep); 
   
%%HERE STARTS EQUATION ASSEMBLING - THIS IS THE BEGINNING OF SETTING UP THE
%%EQUATIONS
   %THE MAIN VARIABLESARE P,F, Zi, Sw So MAKE ADI VARIABLES
   
   fluid=state.fluid
   p=state.pressure;
   F=state.F;
   Zi=state.Zi;
   Sw=state.Sw;
   Xig=state.Xig;
   Xio=state.Xio;
   Ew=state.Ew;
   So=state.So;
   
   cwL=state.cwL
   cwV=state.cwV
   Cw=state.Cw
   
   [p,F,Zi,Sw]=initVariablesADI(p,F,Zi,Sw); %DOES THIS MAKE Xi* ADI ALSO? I BELIEVE IT DOES
   
   fluid0=state0.fluid %THIS ACCOUNTS FOR TEMPERATURE AND PRESSURE, SO I MAY BE BEING SLOPPY?REPETITIVE HERE
   p0=state0.pressure;
   F0=state0.F;
   Zi0=state0.Zi
   Sw0=state0.Sw
   Xig0=state0.Xig;
   Xio0=state0.Xio;
   Ew0=state0.Ew;
   So0=state0.So;
   cwL0=state0.cwL
   cwV0=state0.cwV
   Cw0=state0.Cw %I DONT KNOW IF ALL OF THIS IS NECESSARY
   
[krL,krG]=quadraticRelPerm(So);
bd=bc.dirichlet;
[bc_krL, bc_krG] = quadraticRelperm(bd.So);

g  = norm(gravity);
dz=s.dz; %STILL NEED TO DO THIS IN SETUPSYSTEM

%COMPUTE THE MOBILITIES
mobL=krL./fluid.muL;
mobG=krG./fluid.muG;
bc_mobL   = bc_krL./fluid.muL;
bc_mobG   = bc_krG./fluid.muG;

%COMPUTE UPSTREAM DIRECTION FOR EACH COMPONENT (only including non-water
%components because bravo-dome does, so might need change). Also, I am only
%establishing everything as cell arrays becuase bravo dome does and i think
%they might need to be that way for solvefi.. might change also
dpC = cell(nComp_C, 1);
upC = cell(nComp_C, 1);
for ic = 1:nComp_C
    dpC{ic} = s.p_grad(p) - g*(fluid.components.MW(ic).*dz); %STILL NEED TO SET GRAD STUFF UP IN SETUPSYSTEM AREA
    upC{ic} = (double(dpC{ic})>=0);
end
    dpW = s.p_grad(p) - g*(fluid.components.MW(4).*dz); %COMPONENT 4 IS WATER< I PLAN ON MAKING A FUNCTION THAT REGISTERS WHICH COMPONENTS ARE WHICH SO WE CAN TYPE THEM IN BY NAME INSTEAD
    upW  = (double(dpW)>=0);

%[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,cubic_time]=GI_flash(fluid,thermo,options); %I AM CONCERNED THAT THESE MIGHT NOT COME OUT AS ADI VARIABLES
%I ALREADY KNOW THESE FOR THE INTITIAL FLUID SO I DONT REALLY WANT THIS, BUT IF IT STAYS< I NEED TO TRANSLATE THE TERMS

fluxC=cell(nComp,1); %AGAIN, ONLY CELL BECAUSE BRAVO DOME DOES THAT WAY

    for ic = 1 : nComp
       %%
       % The function |s.faceConcentrations| computes concentration on the faces given
       % cell concentrations, using upwind directions.
       %RESIDUAL FOR NON WATER COMPONENTS
       bc_val = bd.Xig{ic}.*bc_mobG + bd.Xio{ic}.*bc_mobL; 
       fluxC{ic} = s.faceConcentrations(upC{ic}, Xig{ic}.*mobG + Xio{ic}.*mobL, bc_val); %NEED TO ADD SETUP FACE CONCENTRAIONS IN SETUP SYSTEM AREA
       eqs{ic} = (rock.pv/dt).*(F*Zi(ic)-F0*Zi(ic))+ s.div(fluxC{ic}.*rock.T.*dpC{ic});
    end

    % Compute the residual of the mass conservation equation for water.
    bc_val = bd.cwV.*bc_mobG + bd.cwL.*bc_mobL;
    fluxW = s.faceConcentrations(upW, cwV.*mobG + cwG.*mobL, bc_val);%NEED TO ADD IN SETUPCONTROLS
    eqs{nComp + 1} = (rock.pv/dt).*(Ew*Sw - Ew0*Sw0) + s.div(fluxW.*rock.T.*dpW);
    %DONE COMPUTING RESIDUAL FOR WATER
    
    %COMPUTE THE GLOBAL FLOW RESIDUAL
    eqs{nComp+2}=(rock.pv/dt).*(F-F0)+s.div(fluxC{ic}.*rock.T.*dpC{ic}); %THE SECOND TERM IS SAME AS FOR INDIVIDUAL COMPONENTS. THIS MUST CHANGE
    %DONE COMPUTING GLOBAL FLOW EQ
    
    %COMPUTE THE SATURATION RESIDUAL EQUATION
    eqs{nComp+3}=(F-F0)*((So-So0)/(Eo-Eo0)+(Sg-Sg0)/(Eg-Eg0))+Sw-1;
    %DONE COMPUTING THE RESIDUAL FOR SATURATION
    
    %ADD INPUT FLUX
    for ic = 1 : nComp
       eqs{ic}(bc.influx_cells) = eqs{ic}(bc.influx_cells) - bc.C_influx{ic};
    end
    eqs{nComp + 1}(bc.influx_cells) = eqs{nComp + 1}(bc.influx_cells) - bc.water_influx;
    %DONE ADDING INPUT FLUXES
     %THIS IS THE END OF SETTING UP THE EQUATIONS!!!!!!!!!!!!!
   %%NOW CONTINUING WITH SOLVE FI %%I MIGHT WANT TO MAKE THIS EQUATION
   %%ASSEMLE PART SEPERATE. IM GOING TO BED, ILL ASK XIAOMENG TOMORROW
   %%AFTER MEETING WITH THE EXTERNSHIP GROUP
   
   
%}





