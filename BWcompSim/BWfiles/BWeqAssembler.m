%%%HERE STARTS EQUATION ASSEMBLING - THIS IS THE BEGINNING OF SETTING UP THE
%%EQUATIONS
   %THE MAIN VARIABLESARE P,F, Zi, Sw So MAKE ADI VARIABLES
   %
   
   function eqs=BWeqAssembler(tstep, iteration, state0,...
       state,bc,system, ops);
   
   components=system.components;
   nComp=system.nComp;
   rock=system.rock;
   dt=system.options.dt;
   dz=ops.dz;
   p_grad=ops.p_grad;
   div=ops.div;
   faceConcentrations=ops.faceConcentrations;
   

%%   
       opt = struct('Verbose',     mrstVerbose,...
                 'scaling',     [],...
                 'history',     [],  ...
                 'iteration',   -1,  ...
                 'stepOptions', []);
   
  % fluid=state.totalFluid;
   muL=1e-3;
   muG=1e-5;
   muW=9e-4;
   %p_ref =system.p_ref;      
   p_ref=1*atm; %Having problems with 'system' its deleting itself somewhere before solveFi
 %%  
   if tstep~=1 & meta.iteration~=1
  
   Ew=state.Ew
   
   p=state.p;
   F=state.F;
   Zi=state.Zi;
   Sw=state.Sw

[p, F, Zi{:}, Sw]=initVariablesADI(state.p, state.F, state.Zi{:}, state.Sw);
   
   else
      
   Ew=state.Ew;
   p=state.p;
   F=state.F;
   Zi=state.Zi;
   Sw=state.Sw
  

   F0=state0.F;
   Zi0=state0.Zi;
   %Zi0=num2cell(Zi0,1); %THIS IS JUST TO MAKE Zi a 1x3 cell array LIKE C in BravoDome
   Sw0=state0.Sw
   Ew0=state0.Ew;
  
 %%    
[krL,krG]=BWquadraticRelPerm(state.So);
bd=bc.dirichlet;
[bc_krL, bc_krG] = BWquadraticRelPerm(bd.So);
krW=0.5; %TEMPORARY
bc_krW=0.5;%TEMPORARY

g  = norm(gravity);
%dz is already known 

%COMPUTE THE MOBILITIES
mobL=krL./muL;
mobG=krG./muG; %VISC DEFINED IN 
mobW=krW./muW;

rhoW = 1*kilogram/litre;


bc_mobL   = bc_krL./muL;
bc_mobG   = bc_krG./muG;
bc_mobW=bc_krW./muW;

%%
%COMPUTE UPSTREAM DIRECTION FOR EACH COMPONENT (only including non-water
%components because bravo-dome does, so might need change). Also, I am only
%establishing everything as cell arrays becuase bravo dome does and i think
%they might need to be that way for solvefi.. might change also
MW=vertcat(component.MW);
dpC = cell(1,2); %2 is the number of phases and this will be changed/fixed
upC = cell(1,2);
dpC{1}=p_grad(p) - g*(rhoL.*dz);
dpC{2}=p_grad(p) - g*(rhoG.*dz);
for phase = 1:2
    upC{phase} = (double(dpC{phase})>=0);
end
    dpC=dpC'; upC=upC';
    dpW = p_grad(p) - g*(rhoW.*dz); %COMPONENT 4 IS WATER< I PLAN ON MAKING A FUNCTION THAT REGISTERS WHICH COMPONENTS ARE WHICH SO WE CAN TYPE THEM IN BY NAME INSTEAD
    upW  = (double(dpW)>=0);


fluxC=cell(nComp_C,1); %AGAIN, ONLY CELL BECAUSE BRAVO DOME DOES THAT WAY

       %%
       %COMPUTE COMPONENT FLOW RESIDUAL
    for ic = 1 : nComp_C
       bc_val = bd.Xig(ic).*bc_mobG + bd.Xio(ic).*bc_mobL; 
       fluxC{ic} = faceConcentrations(upC{ic}, state.Xig{ic}.*mobG.*state.Eg + state.Xio{ic}.*mobL.*state.Eo, bc_val); %THESE Xi VALUES ARE FOR CELL 1 AND NEED TO BE FIXED
       eqs{ic} = (rock.pv/dt).*(F.*Zi{ic}-F0.*Zi0{ic})+ div(fluxC{ic}.*rock.T.*dpC{ic});
    end
    
    %%
    % Compute the residual of the mass conservation equation for water.
    %NEED TO DEFINE mobW and Bc_mobW
    bc_val = bd.Ew.*bc_mobW;
    fluxW = faceConcentrations(upW, state.Ew.*mobW, bc_val);%NEED TO ADD IN SETUPCONTROLS
    eqs{nComp_C + 1} = (rock.pv/dt).*(Ew.*Sw - Ew0.*Sw0) + div(fluxW.*rock.T.*dpW);

    
    %%
    %COMPUTE THE GLOBAL FLOW RESIDUAL
    bc_val = bd.Eg.*bc_mobG + bd.Eo.*bc_mobL; 
    fluxT = faceConcentrations(upC_total, state.Eg.*mobG + state.Eo.*mobL, bc_val);
    eqs{nComp_C+2}=(rock.pv/dt).*(F-F0)+div(fluxT.*rock.T.*dpC_total); %THE SECOND TERM IS SAME AS FOR INDIVIDUAL COMPONENTS. THIS MUST CHANGE
    %DONE COMPUTING GLOBAL FLOW EQ
    
    %%
    %COMPUTE THE SATURATION RESIDUAL EQUATION
    eqs{nComp_C+3}=(F)*((1-state.V)/(state.Eo)+(state.V)/(state.Eg))+(Sw)-1;
    %DONE COMPUTING THE RESIDUAL FOR SATURATION

        
    %%
    %ADD INPUT FLUX
    for ic = 1 : nComp_C
       eqs{ic}(bc.in.influx_cells) = eqs{ic}(bc.in.influx_cells) - bc.in.C_influx(ic);
    end
    eqs{nComp_C + 1}(bc.in.influx_cells) = eqs{nComp_C + 1}(bc.in.influx_cells) - bc.in.water_influx;
    eqs{nComp_C + 2}(bc.in.influx_cells) = eqs{nComp_C + 2}(bc.in.influx_cells) - bc.in.T_influx;
   
   end
   
    %DONE ADDING INPUT FLUXES
   %NEED TO ADD INPUT FLUXES FOR TOTAL FLOW!!!!!
   
   %IN GENERAL, NEED TO MAKE SURE TOTAL FLOW IS CORRECT
   %Need rho(phase)
   





