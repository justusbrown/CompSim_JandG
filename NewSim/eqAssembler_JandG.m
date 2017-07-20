%%HERE STARTS EQUATION ASSEMBLING - THIS IS THE BEGINNING OF SETTING UP THE
%%EQUATIONS
   %THE MAIN VARIABLESARE P,F, Zi, Sw So MAKE ADI VARIABLES
   
   
   function eqs=eqAssembler_JandG(rock, state0,state,dt,bc,dz,p_grad,div,faceConcentrations);
   %NOT SURE IF ILL KEEP varagin gr 07/20
   
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
   Ew=state.Ew;
   So=state.So;
   
   cwL=state.cwL%NED TO THINK ABOUT ALL OF THESE WATER TERMS gr 07/20
   cwV=state.cwV
   Cw=state.Cw
   
   [p,F,Zi,Sw]=initVariablesADI(p,F,Zi,Sw); %DOES THIS MAKE Xi* ADI ALSO? I BELIEVE IT DOES
   
   fluid0=state0.fluid;%NEEDS THOUGHT: %THIS ACCOUNTS FOR TEMPERATURE AND PRESSURE, SO I MAY BE BEING SLOPPY?REPETITIVE HERE
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
dpC = cell(nComp_C, 1);
upC = cell(nComp_C, 1);
for ic = 1:nComp_C
    dpC{ic} = s.p_grad(p) - g*(fluid.components.MW(ic).*dz); %STILL NEED TO SET GRAD STUFF UP IN SETUPSYSTEM AREA
    upC{ic} = (double(dpC{ic})>=0);
end
    dpW = p_grad(p) - g*(fluid.components.MW(4).*dz); %COMPONENT 4 IS WATER< I PLAN ON MAKING A FUNCTION THAT REGISTERS WHICH COMPONENTS ARE WHICH SO WE CAN TYPE THEM IN BY NAME INSTEAD
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
       fluxC{ic} = faceConcentrations(upC{ic}, Xig{ic}.*mobG + Xio{ic}.*mobL, bc_val);
       eqs{ic} = (rock.pv/dt).*(F*Zi(ic)-F0*Zi0(ic))+ div(fluxC{ic}.*rock.T.*dpC{ic});
    end

    % Compute the residual of the mass conservation equation for water.
    bc_val = bd.cwV.*bc_mobG + bd.cwL.*bc_mobL;
    fluxW = faceConcentrations(upW, cwV.*mobG + cwG.*mobL, bc_val);%NEED TO ADD IN SETUPCONTROLS
    eqs{nComp + 1} = (rock.pv/dt).*(Ew*Sw - Ew0*Sw0) + div(fluxW.*rock.T.*dpW);
    %DONE COMPUTING RESIDUAL FOR WATER
    %%STILL NEED Ew gr 07/20
    %I GUESS WE DO NEED cw stuff gr 07/20
    
    %COMPUTE THE GLOBAL FLOW RESIDUAL
           bc_val = bd.Xig{ic}.*bc_mobG + bd.Xio{ic}.*bc_mobL; 
       fluxT = faceConcentrations(upC{ic}, Xig{ic}.*mobG + Xio{ic}.*mobL, bc_val);
       %%THIS NEEDS THOUGHT gr 07/20
    eqs{nComp+2}=(rock.pv/dt).*(F-F0)+div(fluxT{ic}.*rock.T.*dpC{ic}); %THE SECOND TERM IS SAME AS FOR INDIVIDUAL COMPONENTS. THIS MUST CHANGE
    %DONE COMPUTING GLOBAL FLOW EQ
    
    %COMPUTE THE SATURATION RESIDUAL EQUATION
    eqs{nComp+3}=(F-F0)*((So-So0)/(Eo-Eo0)+(Sg-Sg0)/(Eg-Eg0))+(Sw-Sw0)-1;
    %DONE COMPUTING THE RESIDUAL FOR SATURATION
    %NEED TO KNOW WHAT IS SUBTRACTED AND WHAT ISNT. PROBABLY PRETTY EASY TO
    %LOOK UP gr 07/20
    
    
    %DO JUST PRIMARY VARIABLES NEED TO BE (X-X0)???? gr 07/20
    
    %ADD INPUT FLUX
    for ic = 1 : nComp
       eqs{ic}(bc.influx_cells) = eqs{ic}(bc.influx_cells) - bc.C_influx{ic};
    end
    eqs{nComp + 1}(bc.influx_cells) = eqs{nComp + 1}(bc.influx_cells) - bc.water_influx;
    
   end
    %DONE ADDING INPUT FLUXES
     %THIS IS THE END OF SETTING UP THE EQUATIONS!!!!!!!!!!!!!
   %%NOW CONTINUING WITH SOLVE FI %%I MIGHT WANT TO MAKE THIS EQUATION
   %%ASSEMLE PART SEPERATE. IM GOING TO BED, ILL ASK XIAOMENG TOMORROW
   %%AFTER MEETING WITH THE EXTERNSHIP GROUP
   
   





