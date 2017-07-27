%%%HERE STARTS EQUATION ASSEMBLING - THIS IS THE BEGINNING OF SETTING UP THE
%%EQUATIONS
   %THE MAIN VARIABLESARE P,F, Zi, Sw So MAKE ADI VARIABLES
   %
   
   function eqs=eqAssembler_JandG(tstep, iteration, component, rock, state0,state,dt,bc,dz,p_grad,div,faceConcentrations);
   %should nComp_C just be an input?
   nComp_C=3; %jb 7/24
   
       opt = struct('Verbose',     mrstVerbose,...
                 'scaling',     [],...
                 'history',     [],  ...
                 'iteration',   -1,  ...
                 'stepOptions', []);

   %THESE ARE ALL INSIDE STATE... LIFE MAY BE EASIER IF WE DONT REDEFINE HERE 
    
   %{
   p=state.totalFluid.pressure;
   F=state.totalFluid.F;
   Zi=state.totalFluid.Zi;
   %Zi=num2cell(Zi,1); %THIS IS JUST TO MAKE Zi a 1x3 cell array LIKE C in BravoDome
   Sw=state.totalFluid.Sw;
   Xig=state.totalFluid.Xig;
   Xio=state.totalFluid.Xio;
   Xwv=state.totalFluid.Xwv;
   Xwl=state.totalFluid.Xwl;
   Eo=state.totalFluid.Eo;
   Eg=state.totalFluid.Eg;
   Ew=state.totalFluid.Ew;
   So=state.totalFluid.So;
   V=state.totalFluid.V;

   %}
   
   fluid=state.totalFluid;
   muL=1e-3;
   muG=1e-5;
   p_ref = 1*atm;       % Reference pressure
   
   %combined=[fluid(1:rock.G.cells.num)];
   
   %[Xig,Xio,Xwv,Xwl,Eo,Eg,Ew,So,V]=variableCall_JandG(rock,state);
   %Xig=Xig';Xio=Xio';Xwv=Xwv';Xwl=Xwl';Eo=Eo';Eg=Eg';Ew=Ew';So=So';V=V';
   
   %[p,F,Sw,Zi]=primeVars_JandG(rock, state);
   %p=p'; F=F'; Sw=Sw'; Zi=Zi';
   
   %I am only explicitly defining things in finite differences or primaries
  if tstep~=1 & meta.iteration~=1
   Ew=state.Ew
   %PRIME VARS
   p=state.p;
   F=state.F;
   Zi=state.Zi;
   Sw=state.Sw

[p, F, Zi{:}, Sw]=initVariablesADI(state.p, state.F, state.Zi{:}, state.Sw);
   
  else
   Ew=state.Ew;
   %PRIME VARS
   p=state.p;
   F=state.F;
   Zi=state.Zi;
   Sw=state.Sw
  

   %DOES ALL OF THIS NEED TO BE INCLUDED, COMMENTED OUT WHAT ISN'T USED IN
   %THIS SHEET. JB 7/21
   %fluid0=state0.fluid;%NEEDS THOUGHT: %THIS ACCOUNTS FOR TEMPERATURE AND PRESSURE, SO I MAY BE BEING SLOPPY?REPETITIVE HERE
   F0=state0.F;%USED LATER
   Zi0=state0.Zi;%USED LATER
   %Zi0=num2cell(Zi0,1); %THIS IS JUST TO MAKE Zi a 1x3 cell array LIKE C in BravoDome
   Sw0=state0.Sw%USED LATER
   Ew0=state0.Ew;%USED LATER
   %p0=state0.pressure;
   %So0=state0.So;
   %Xig0=state0.Xig;
   %Xio0=state0.Xio;
   
   %STILL NEED TO BE COMPUTED
   %cwL0=state0.cwL;
   %cwV0=state0.cwV;
   %Cw0=state0.Cw; %I DONT KNOW IF ALL OF THIS IS NECESSARY

   %added this function under the calling file jb 7/24
[krL,krG]=quadraticRelPerm_JandG(state.So);
bd=bc.dirichlet;
[bc_krL, bc_krG] = quadraticRelPerm_JandG(bd.So);

g  = norm(gravity);
%dz is already known 

%COMPUTE THE MOBILITIES
mobL=krL./muL;
mobG=krG./muG; %VISC DEFINED IN 
bc_mobL   = bc_krL./muL;
bc_mobG   = bc_krG./muG;

%COMPUTE UPSTREAM DIRECTION FOR EACH COMPONENT (only including non-water
%components because bravo-dome does, so might need change). Also, I am only
%establishing everything as cell arrays becuase bravo dome does and i think
%they might need to be that way for solvefi.. might change also
MW=vertcat(component.MW);
dpC = cell(1,nComp_C);
upC = cell(1,nComp_C);
for ic = 1:nComp_C
    dpC{ic} = p_grad(p) - g*(MW(ic).*dz); 
    upC{ic} = (double(dpC{ic})>=0);
end
    dpC=dpC'; upC=upC';
    dpW = p_grad(p) - g*(MW(4).*dz); %COMPONENT 4 IS WATER< I PLAN ON MAKING A FUNCTION THAT REGISTERS WHICH COMPONENTS ARE WHICH SO WE CAN TYPE THEM IN BY NAME INSTEAD
    upW  = (double(dpW)>=0);

%[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,cubic_time]=GI_flash(fluid,thermo,options); %I AM CONCERNED THAT THESE MIGHT NOT COME OUT AS ADI VARIABLES
%I ALREADY KNOW THESE FOR THE INTITIAL FLUID SO I DONT REALLY WANT THIS, BUT IF IT STAYS< I NEED TO TRANSLATE THE TERMS

fluxC=cell(nComp_C,1); %AGAIN, ONLY CELL BECAUSE BRAVO DOME DOES THAT WAY

    for ic = 1 : nComp_C
       %%
       % The function |s.faceConcentrations| computes concentration on the faces given
       % cell concentrations, using upwind directions.
       %RESIDUAL FOR NON WATER COMPONENTS
       bc_val = bd.Xig(ic).*bc_mobG + bd.Xio(ic).*bc_mobL; 
       fluxC{ic} = faceConcentrations(upC{ic}, state.Xig{ic}.*mobG + state.Xio{ic}.*mobL, bc_val); %THESE Xi VALUES ARE FOR CELL 1 AND NEED TO BE FIXED
       eqs{ic} = (rock.pv/dt).*(F.*Zi{ic}-F0.*Zi0{ic})+ div(fluxC{ic}.*rock.T.*dpC{ic});
    end

    % Compute the residual of the mass conservation equation for water.
    bc_val = bd.Xwv.*bc_mobG + bd.Xwl.*bc_mobL;
    fluxW = faceConcentrations(upW, state.Xwv.*mobG + state.Xwl.*mobL, bc_val);%NEED TO ADD IN SETUPCONTROLS
    eqs{nComp_C + 1} = (rock.pv/dt).*(Ew.*Sw - Ew0.*Sw0) + div(fluxW.*rock.T.*dpW);
    %DONE COMPUTING RESIDUAL FOR WATER
    %%STILL NEED Ew gr 07/20
    %I GUESS WE DO NEED cw stuff gr 07/20
    %DO WE HAVE CWV AND SUCH, SEEMS LIKE A LOT OF THIS IS DONE IN
    %COMPUTEWATER IN BRAVO, ANY THOUGHTS? JB 7/21
    
    %%
    %COMPUTE THE GLOBAL FLOW RESIDUAL
    %
    %First need to define things: avg. MW, global dpC and global upC
    %gr 07/20
    %avgMW=sum(fluid.mole_fraction(1:3)'.*MW(1:3));
    %avgMW=getAvgMW(components, Zi);
    dpC_total = p_grad(p) %- g*(avgMW.*dz); 
    upC_total = (double(dpC_total)>=0);

    bc_val = bd.V.*bc_mobG + (1-bd.V).*bc_mobL; 
    fluxT = faceConcentrations(upC_total, state.V.*mobG + (1-state.V).*mobL, bc_val);
       %CHANGED Xia to V an 1-V. dpC{ic} NEEDS THOUGHT HERE AND I need to
       %find out if faceConc can take in a single value for varagin 2
       %%THIS NEEDS THOUGHT gr 07/20
    eqs{nComp_C+2}=(rock.pv/dt).*(F-F0)+div(fluxT.*rock.T.*dpC_total); %THE SECOND TERM IS SAME AS FOR INDIVIDUAL COMPONENTS. THIS MUST CHANGE
    %DONE COMPUTING GLOBAL FLOW EQ
    
    %COMPUTE THE SATURATION RESIDUAL EQUATION
    eqs{nComp_C+3}=(F)*((1-state.V)/(state.Eo)+(state.V)/(state.Eg))+(Sw)-1;
    %DONE COMPUTING THE RESIDUAL FOR SATURATION

        
    %ADD INPUT FLUX
    for ic = 1 : nComp_C
       eqs{ic}(bc.in.influx_cells) = eqs{ic}(bc.in.influx_cells) - bc.in.C_influx(ic);
    end
    eqs{nComp_C + 1}(bc.in.influx_cells) = eqs{nComp_C + 1}(bc.in.influx_cells) - bc.in.water_influx;
    
   end
    %DONE ADDING INPUT FLUXES
   %NEED TO ADD INPUT FLUXES FOR TOTAL FLOW!!!!!
   
   %IN GENERAL, NEED TO MAKE SURE TOTAL FLOW IS CORRECT
   





