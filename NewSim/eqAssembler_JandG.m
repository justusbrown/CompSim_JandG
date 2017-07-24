%%%HERE STARTS EQUATION ASSEMBLING - THIS IS THE BEGINNING OF SETTING UP THE
%%EQUATIONS
   %THE MAIN VARIABLESARE P,F, Zi, Sw So MAKE ADI VARIABLES
   %
   
   function eqs=eqAssembler_JandG( rock, state0,state,dt,bc,dz,p_grad,div,faceConcentrations);
   %should nComp_C just be an input?
   nComp_C=3; %jb 7/24
   
       opt = struct('Verbose',     mrstVerbose,...
                 'scaling',     [],...
                 'history',     [],  ...
                 'iteration',   -1,  ...
                 'stepOptions', []);

    
    
   fluid=state.fluid; %NEED TO THINK ABOUT THIS gr 07/20
   p=state.pressure;
   F=state.F;
   Zi=state.Zi;
   Sw=state.Sw;
   Xig=state.Xig;
   Xio=state.Xio;
   Xwv=state.Xwv;
   Xwl=state.Xwl;
   Ew=state.Ew;
   So=state.So;

   
   [p,F,Zi,Sw]=initVariablesADI(p,F,Zi,Sw); %DOES THIS MAKE Xi* ADI ALSO? I BELIEVE IT DOES
   
   %DOES ALL OF THIS NEED TO BE INCLUDED, COMMENTED OUT WHAT ISN'T USED IN
   %THIS SHEET. JB 7/21
   fluid0=state0.fluid;%NEEDS THOUGHT: %THIS ACCOUNTS FOR TEMPERATURE AND PRESSURE, SO I MAY BE BEING SLOPPY?REPETITIVE HERE
   F0=state0.F;%USED LATER
   Zi0=state0.Zi%USED LATER
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
[krL,krG]=quadraticRelPerm_JandG(So);
bd=bc.dirichlet;
[bc_krL, bc_krG] = quadraticRelPerm_JandG(bd.So);

g  = norm(gravity);
%dz is already known 

%COMPUTE THE MOBILITIES
mobL=krL./fluid.muL;
mobG=krG./fluid.muG; %VISC DEFINED IN 
bc_mobL   = bc_krL./fluid.muL;
bc_mobG   = bc_krG./fluid.muG;

%COMPUTE UPSTREAM DIRECTION FOR EACH COMPONENT (only including non-water
%components because bravo-dome does, so might need change). Also, I am only
%establishing everything as cell arrays becuase bravo dome does and i think
%they might need to be that way for solvefi.. might change also
MW=vertcat(fluid.components.MW);
dpC = cell(1,nComp_C);
upC = cell(1,nComp_C);
for ic = 1:nComp_C
    dpC{ic} = p_grad(p) - g*(MW(ic).*dz); %STILL NEED TO SET GRAD STUFF UP IN SETUPSYSTEM AREA
    upC{ic} = (double(dpC{ic})>=0);
end
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
       bc_val = bd.Xig{ic}.*bc_mobG + bd.Xio{ic}.*bc_mobL; 
       fluxC{ic} = faceConcentrations(upC{ic}, Xig{ic}.*mobG + Xio{ic}.*mobL, bc_val);
       eqs{ic} = (rock.pv/dt).*(F*Zi(ic)-F0*Zi0(ic))+ div(fluxC{ic}.*rock.T.*dpC{ic});
    end

    % Compute the residual of the mass conservation equation for water.
    bc_val = bd.Xwv.*bc_mobG + bd.Xwl.*bc_mobL;
    fluxW = faceConcentrations(upW, Xwv.*mobG + Xwl.*mobL, bc_val);%NEED TO ADD IN SETUPCONTROLS
    eqs{nComp + 1} = (rock.pv/dt).*(Ew*Sw - Ew0*Sw0) + div(fluxW.*rock.T.*dpW);
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
    avgMW=sum(fluid.mole_fraction(1:3).*fluid.component.MW(1:3));
    dpC_total = s.p_grad(p) - g*(avgMW.*dz); 
    upC_total = (double(dpC_total)>=0);

    bc_val = bd.V.*bc_mobG + (1-bd.V).*bc_mobL; 
    fluxT = faceConcentrations(upC_total, V*mobG + (1-V).*mobL, bc_val);
       %CHANGED Xia to V an 1-V. dpC{ic} NEEDS THOUGHT HERE AND I need to
       %find out if faceConc can take in a single value for varagin 2
       %%THIS NEEDS THOUGHT gr 07/20
    eqs{nComp+2}=(rock.pv/dt).*(F-F0)+div(fluxT.*rock.T.*dpC_total); %THE SECOND TERM IS SAME AS FOR INDIVIDUAL COMPONENTS. THIS MUST CHANGE
    %DONE COMPUTING GLOBAL FLOW EQ
    
    %COMPUTE THE SATURATION RESIDUAL EQUATION
    eqs{nComp+3}=(F)*((1-V)/(Eo)+(V)/(Eg))+(Sw)-1;
    %DONE COMPUTING THE RESIDUAL FOR SATURATION

        
    %ADD INPUT FLUX
    for ic = 1 : nComp
       eqs{ic}(bc.influx_cells) = eqs{ic}(bc.influx_cells) - bc.C_influx{ic};
    end
    eqs{nComp + 1}(bc.influx_cells) = eqs{nComp + 1}(bc.influx_cells) - bc.water_influx;
    
   end
    %DONE ADDING INPUT FLUXES
   
   





