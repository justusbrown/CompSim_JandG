function eqs = BWeqAssembler(state0, state, bc, system, ops)

 opt = struct('Verbose',     mrstVerbose,...
                 'scaling',     [],...
                 'history',     [],  ...
                 'iteration',   -1,  ...
                 'stepOptions', []);   

   components=system.components;
   nComp=system.nComp;
   nPhase=2;
   rock=system.rock;
   dt=system.options.dt;
   dz=ops.dz;
   p_grad=ops.p_grad;
   div=ops.div;
   faceConcentrations=ops.faceConcentrations;

%Clear/initialize the eqs variables
  eqs=cell(1,nComp+2);
  
 %Viscosities will not be here in final
   muL=1e-3;
   muG=1e-5;
   muW=9e-4;
   
 %ADI/Primary Variables
 p=state.p;
 m_i=state.m_i;
 m_w=state.m_w;
 
 [p, m_i{:}, m_w]=initVariablesADI(state.p, state.m_i{:}, state.m_w);
 
 %ADDED THIS!!!!! JB 8/15
 for i=1:system.nCell
 m_iSum(i) = m_i{1}(i) + m_i{2}(i) + m_i{3}(i) + m_i{4}(i) + m_i{5}(i) + m_i{6}(i) ; %assuming constant for each cell
 end
 m_iSum=m_iSum(:)';
 
 %Finite Difference Variables
 m_i0=state0.m_i;
 m_w0=state0.m_w;
 
  %% Setup the flux terms   %HOW DID YOU GET THIS BOLDED?
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

MW=vertcat(components.MW);
rhoL=state.rhoL;
rhoG=state.rhoG;

fz=rock.G.faces.centroids(:,3);
drhoL=ops.grad(rhoL, fz);
drhoG=ops.grad(rhoG, fz);


dpC = cell(1,2);
upC = cell(1,2);
dpC{1}=p_grad(p) - g*(drhoL.*dz);
dpC{2}=p_grad(p) - g*(drhoG.*dz);
for phase = 1:nPhase
    upC{phase} = (double(dpC{phase})>=0);
end
    dpC=dpC'; upC=upC';
    dpW = p_grad(p) - g*(rhoW.*dz); 
    upW  = (double(dpW)>=0);
 
fluxL=cell(nComp,1);
fluxG=cell(nComp,1);
fluxC=cell(nComp,1);

%%Compute the component flow residual
for ic = 1 : nComp
   bc_valG = bd.Xig(ic).*bc_mobG.*bd.Eg;
   bc_valL= bd.Xio(ic).*bc_mobL.*bd.Eo; 
   valL=state.Xio{ic}.*mobL.*state.Eo;
   valG=state.Xig{ic}.*mobG.*state.Eg;
   fluxL{ic} = faceConcentrations(upC{1}, valL, bc_valL);% CHANGED TO BWFACECONCENTRATIONS
   fluxG{ic}= faceConcentrations(upC{2}, valG, bc_valG);
   fluxC{ic}=fluxL{ic}.*dpC{1}+fluxG{ic}.*dpC{2};
   eqs{ic} = (rock.pv/dt).*(m_i{ic}-m_i0{ic})+ div(fluxC{ic}.*rock.T);
end

%%Compute the water flow residual
bc_val = bd.Ew.*bc_mobW;
fluxW = faceConcentrations(upW, state.Ew.*mobW, bc_val);%NEED TO ADD IN SETUPCONTROLS
eqs{nComp + 1} = (rock.pv/dt).*(m_w - m_w0) + div(fluxW.*rock.T.*dpW);

%%Compute the saturation residual 
%THE WAY IT WAS:eqs{nComp+3}=state.m_sum*(1-state.V)/state.Eo + state.m_sum*state.V/state.Eg + m_w/state.Ew-1;
%CHANGED THIS, MAY BE COMPLETE GARBAGE JB 8/15
eqs{nComp+2}=m_iSum.*(1-state.V)./state.Eo + m_iSum.*state.V./state.Eg + m_w./state.Ew-1;

%%Add the input fluxes
for ic = 1 : nComp
       eqs{ic}(bc.in.influx_cells) = eqs{ic}(bc.in.influx_cells) - bc.in.C_influx(ic);
end
    eqs{nComp + 1}(bc.in.influx_cells) = eqs{nComp + 1}(bc.in.influx_cells) - bc.in.water_influx;

   
end



 