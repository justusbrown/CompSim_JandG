%%
%This is where the user where input their data to run our simulator.
%
%USER INPUTS
%
%***All of the objects/data that are setup in this function call are listed in
%detail at the end of the function

%FUNCTION inputData

function [system]=BWinputData()

system=struct();
%
%%
%%%ENTER THE ROCK DATA (GOOD)
%
%%Enter the desired 3x3 grid dimmensions 
x=293.3*(0:1:9)*ft;
y=293.3*(0:1:9)*ft;
z1=30*(0:1:2);
z2=z1(3)+50*(1:1:2);
z=horzcat(z1,z2)*ft;
G = tensorGrid(x,y,z);
G = computeGeometry(G);
%%Enter the permeability
rock.perm=repmat(100*milli*darcy, [G.cells.num, 1]);
%%Enter the porosity
rock.poro = repmat(0.13, [G.cells.num,1]);
rock.pv=poreVolume(G,rock);
rock.T=computeTrans(G,rock);
rock.G=G
system.Cells=rock.G.cells;
system.nCell=rock.G.cells.num;


%%
%%%ENTER THE PVTi OPTIONS (GOOD)
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
%%%!GOOD
%%%AND CHOICE OF EOS/MIXING RULES (included in addThermo for our sake)

%Enter the influx pressure in Pa
system.influx_p=10e6;
%Enter the outflux pressure in Pa
system.outflux_p=8e6;
%Assuming system pressure=outflux pressure
for i=1:system.nCell
    system.Cells.pressure(i)=system.outflux_p;
end
%Enter the influx  in m^3/s
system.influx_rate = 1000/day  
%Enter the temperature in Kelvin
system.Temp=(200+460)*Rankine;
system.thermo=addThermo();
system.thermo.EOS=@PREOS;
system.thermo.vp_water=BWvaporPressure(system.Temp);
system.thermo.R=8.3145;

%%
%ENTER FLUID COMPONENTS !MAYBE CHANGE
namecomponents=['CO2', 'N2', 'C1', 'C2', 'C3', 'C4-6', 'C7+1','C7+2', 'C7+3'];
system.nComp=length(namecomponents);

%%

%CRITICAL PROPERTIES !GOOD
components.Pc=[1071.33111,492.31265,667.78170,708.34238,618.69739,514.92549,410.74956,247.56341,160.41589]*psia;
components.Tc=[548.46000,227.16000,343.08000,549.77400,665.64000,806.54054,838.11282,1058.03863,1291.89071]*Rankine;
components.Zc=[0.27408,0.29115,0.28473,0.28463,0.27748,0.2764,0.2612,0.22706,0.20137];
%%
%COMPONENT PROPERTIES
components.acentric_factor=[0.225,0.04,0.013,0.0986,0.1524,0.21575,0.3123,0.5567,0.91692];
components.MW=[16.04,44.10,86.18,149.29,206.00,282.00]*milli; %FIX
components.OMEGA_A=[0.4572355,0.4572355,0.534021,0.4572355,0.4572355,0.4572355,0.6373344,0.6373344,0.6373344];
components.OMEGA_B=[0.0777961,0.0777961,0.0777961,0.0777961,0.0777961,0.0777961,0.0872878,0.0872878,0.0872878];
components.refRHO=[48.50653,50.19209,26.53189,34.21053,36.33308,37.87047,45.60035,50.88507,55.89861]*lbm_ft3;
%%
%SATURATION DATA NEW

%WATER
SAT_W=[0.16,0.18,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.760,0.8,0.84,0.88,0.92,0.96,1.00];
Kr_W=[0,0,0.002,0.010,0.020,0.033,0.049,0.066,0.090,0.119,0.150,0.186,0.227,0.277,0.330,0.390,0.462,0.540,0.620,0.710,0.800,0.900,1.000];
CapP_W=[50,41,32,21,15.5,12.0,9.2,7.0,5.3,4.2,3.4,2.7,2.1,1.7,1.3,1.0,0.7,0.5,0.4,0.3,0.2,0.1,0.0];

%GAS
Sat_G=[0.00,0.04,0.08,0.12,0.16,0.18,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.760,0.8,0.84];
Kr_G=[0.000,0.005,0.013,0.026,0.040,0.058,0.078,0.100,0.126,0.156,0.187,0.222,0.260,0.300,0.349,0.400,0.450,0.505,0.562,0.620,0.680,0.740];
CapP_G=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1];

%OIL
Sat_O=[0.00,0.04,0.08,0.12,0.16,0.18,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.760,0.8,0.84];
Kr_O=[0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.005,0.012,0.024,0.040,0.060,0.082,0.112,0.150,0.196,0.315,0.400,0.513,0.650,0.800];
CapP_O=[0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.005,0.012,0.024,0.040,0.060,0.082,0.112,0.150,0.196,0.315,0.400,0.513,0.650,0.800];
%%
%ENTER INITIAL MOLE FRACTIONS

Zi=[.01210,.01940,.65990 ,.08690,.05910, .09670,.04745,.01515,.00330];
m_i=[.01210,.01940,.65990 ,.08690,.05910, .09670,.04745,.01515,.00330];%need to do something here

%Assign the initialFluid to the entire system using addBWfluid
system.components=components;
totalFluid=addBWfluid(system, Zi, m_i);
%Initialize the influx, outflux fluida
influxFluid=addMixture(system.components,system.Temp,system.influx_p);
outfluxFluid=addMixture(system.components,system.Temp,system.outflux_p);
%Enter the Influx Fluid's mole fraction
influxFluid.Zi=[0.4,0.03,0.17,0.1,0.25,0.05];%NEEDS TO CHANGE JB 9/1
%Enter the Outflux Fluid's mole fraction
outfluxFluid.Zi=[0.5,0.03,0.07,0.2,0.15,0.05];%NEEDS TO CHANGE JB 9/1
%%
%%%Enter the nonlinear solver parameters and ***cellwise***
%
options.maxIterations=50;
options.nonlinear=setBWnonlinearSolverParameters(options.maxIterations);
options.cellwise=1:5;

%%
%%%Enter the time options for the solver
%
%Enter the time step
options.dt = 200*day;    
%Enter the total time
options.total_time = 10*365*day;  %CHANGED
options.steps      = options.dt*ones(floor(options.total_time/options.dt), 1); 
options.t          = cumsum(options.steps); 

%%
%GROUP EVERYTHING INTO OVERARCHING SYSTEM
system.options=options;
system.totalFluid=totalFluid;
system.influxFluid=influxFluid;
system.outfluxFluid=outfluxFluid;
system.rock=rock;
system.cl=4.4e-5/atm; % Compressibility
system.p_ref = 1*atm;      % Reference pressure
%water info
mmH  = 1.00794*gram;  % molar mass of Hydrogen
mmO  = 15.9994*gram;  % molar mass of Oxygen
mmW = 2*mmH + mmO;  % molar mass of H20
system.litre = 1e-3*meter^3;
rho = 1*kilogram/system.litre;
system.mv_water=mmW/rho;



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

