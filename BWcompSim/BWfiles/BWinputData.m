%%
%This is where the user where input their data to run our simulator.
%
%USER INPUTS
%
%***All of the objects/data that are setup in this function call are listed in
%detail at the end of the function

%FUNCTION inputData

function [G,rock,options,thermo,influxFluid,outfluxFluid,initialFluid, influx_rate, system]=inputData()


%
%%
%%%ENTER THE ROCK DATA
%
%%Enter the desired 3x3 grid dimmensions
G = cartGrid([10, 10, 4]); %CHANGED!!!!
G = computeGeometry(G);
%%Enter the permeability
rock.perm=repmat(50*milli*darcy, [G.cells.num, 1]);
%%Enter the porosity
rock.poro = repmat(0.1, [G.cells.num,1]);
rock.pv=poreVolume(G,rock);
rock.T=computeTrans(G,rock);
rock.G=G

%%
%%%ENTER THE PVTi OPTIONS
%
%Enter the convergance tolerance for fugacity
options.convergence_eps = 1e-12;
%Enter the trivial shift for bisection algorithm
options.trivial_eps = 1e-3;
%Enter the max number of RR iterations using Newton's Methods
options.RRiteration = 200;
%Enter the max number of fugacity updates
options.max_outer_loop = 1000;

%%
%%%ENTER THE INFLUX AND OUTFLUX PRESSURE DATA, INFLUX RATE, TEMPERATURE,
%%%AND CHOICE OF EOS/MIXING RULES (included in addThermo for our sake)
%
%WHERE'S THIS INFO IN THE DECK!!!
%Enter the influx pressure in Pa
influx_p=10e6;
%Enter the outflux pressure in Pa
outflux_p=8e6;
%Enter the influx  in m^3/s
influx_rate = 1000/day  
%Enter the temperature in Kelvin

%KEPT IT AS THEIR UNITS
Temp=((160-32)*(5/9))+273.15;
thermo=addThermo();
thermo.EOS=@PREOS;
thermo.vp_water=vaporPressure(Temp);

%%
%%%ENTER THE COMPONENENTS AND MOLE FRACTIONS FOR INFLUX, OUTFLUX, AND
%%%INITIAL FLUIDS
%%NOTE: Here, fluid properties will be input, but for development sake, our
%%tables will be used
%
%Enter all fluid components NEW!!!!!
[components, comp_flag]=addComponents({'CH4','C3H8','C6H14','C10H22','C15H32','C20H42'});
%Initialize the influx, outflux, and initial fluids
influxFluid=addMixture(components,Temp,influx_p);
outfluxFluid=addMixture(components,Temp,outflux_p);
initialFluid=addMixture(components,Temp,outflux_p);
%Enter the Influx Fluid's mole fraction
%JUST MADE THIS UP NEW!!!!!!a
influxFluid.Zi=[0.4,0.03,0.17,0.1,0.25,0.05];
%Enter the Outflux Fluid's mole fraction
outfluxFluid.Zi=[0.5,0.03,0.07,0.2,0.15,0.05];
%Enter the Initial Fluid's mole fraction
initialFluid.Zi=[0.5,0.03,0.07,0.2,0.15,0.05];
%%NOTE THAT THE MOLE FRACTION ENTERED CORRESPONDS TO THE ORDER OF
%%COMPONENTS ENTERED

%%
%NEW!!!!!!!! COMMENTED OUT BECAUSE WE DON'T HAVE ANYTHING CALLED FLUID
%CRITICAL PROPERTIES
%CHANGE THIS TO WHATEVER STRUCT THIS NEEDS TO BE
%NOT SURE ABOUT THESE UNITS YET
%fluid.Pcrit=[667.8,616.3,436.9,304.0,200.0,162.0]
%fluid.Tcrit=[343.0,665.7,913.4,1011.8,1270.0,1380.0]
%fluid.Zcrit=[.290,.277,.264,.257,.245,.235]
%%
%NEW!!!!!!!
%fluid.AcenF=[.013,.1524,.3007,.4885,.6500,.8500]
%fluid.MMW=[16.04,44.10,86.18,149.29,206.00,282.00]
%%
%%%Enter the nonlinear solver parameters and ***cellwise***
%
maxIterations=50;
nonlinear=setBWnonlinearSolverParameters(maxIterations);
cellwise=1:5;

%%
%%%Enter the time options for the solver
%
%Enter the time step
dt = 200*day;    
%Enter the total time
total_time = 10*365*day;  %CHANGED
steps      = dt*ones(floor(total_time/dt), 1); 
t          = cumsum(steps); 

%%
%GROUP EVERYTHING INTO OVERARCHING SYSTEM
system.R=8.3145;
system.Temp=Temp;
system.vp=thermo.vp_water;
system.fluid=[influxFluid,outfluxFluid,initialFluid];
system.components=components;
system.cl=4.4e-5/atm; % Compressibility
system.p_ref = 1*atm;      % Reference pressure
influxFluid.call=1; %These are setup to avoid confusion when referenceing the fluid vector
outfluxFluid.call=2;
initialFluid.call=3;
system.nComp=numel(components);
%water info
mmH  = 1.00794*gram;  % molar mass of Hydrogen
mmO  = 15.9994*gram;  % molar mass of Oxygen
mmW = 2*mmH + mmO;  % molar mass of H20
litre = 1e-3*meter^3;
rho = 1*kilogram/litre;
system.mv = mmW/rho; %molar volume of water
system.nonlinear=nonlinear;
system.cellwise=cellwise;
system.dt=dt;
system.total_time=total_time;
system.steps=steps;
system.t=t;


end


%%
%***HERE IS A LIST OF EVERYTHING THAT IS ESTABLISHED THROUGH THIS FUNCTION

%%
%{
G  struct with fields:
      cells: [1×1 struct]
       faces: [1×1 struct]
      nodes: [1×1 struct]
      cartDims: [6 6 6]
      type: {'tensorGrid'  'cartGrid'  'computeGeometry'}
      griddim: 3

G.cells  struct with fields:
       num: 216
      facePos: [217×1 double]
     indexMap: [216×1 double]
        faces: [1296×2 double]
      volumes: [216×1 double]
    centroids: [216×3 double]

G.faces struct with fields:
         num: 756
      nodePos: [757×1 double]
    neighbors: [756×2 double]
          tag: [756×1 double]
        nodes: [3024×1 double]
        areas: [756×1 double]
      normals: [756×3 double]
    centroids: [756×3 double]

G.nodes struct with fields:
       num: 343
    coords: [343×3 double]
%}

%{
rock struct with fields:
     perm: [216×1 double]
     poro: [216×1 double]
       pv: [216×1 double]
        T: [756×1 double]
        G: [1×1 struct]
%}

%{
options  struct with fields:
    convergence_eps: 1.0000e-12
        trivial_eps: 1.0000e-03
        RRiteration: 200
     max_outer_loop: 1000
%}

%{
Thermo  struct with fields:
         mixingrule: 1
                EOS: @PREOS
              phase: 1
    fugacity_switch: 1 %%%%vp_water needs to be added here
%}

%{
influxFluid, outfluxFluid, initialFluid struct with fields: 
              bip: [1×1 struct]
       components: [1×4 struct]
    Zi: [0.2000 0.4000 0.3000 0.1000]
         pressure: 10000000
      temperature: 303.1500

influxFluid, outfluxFluid, initialFluid .components 1×4 (nComp) struct array with fields:
   name
    formula
    MW
    Tc
    Pc
    Vc
    Zc
    acentric_factor
    Psat_eq
    Psat_coefs
    PsatTrange
    Psatrange
    dh_vap_eq
    dh_vap_coefs
    dh_vap_Trange
    dh_vap_range
    cp_liq_eq
    cp_liq_coefs
    cp_liq_Trange
    cp_liq_range
    cp_ig_eq
    cp_ig_coefs
    cp_ig_Trange
    cp_ig_range
    dhf_ig
    dgf_ig
    ds_ig
    dh_comb
%}

