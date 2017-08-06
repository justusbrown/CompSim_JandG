function [state, convergence] = BWsolveFI(tstep, system, ops, thermo, rock, state0, bc, equation, options)

    [components, dt, dz, p_grad, div, faceConcentrations]=deal(system.components,...
         system.dt, ops.dz, ops.p_grad, ops.div, ops.faceConcentrations);
   
   meta = struct('converged'  , false                       , ...
                 'stopped'    , false                       , ...
                 'relax'      , system.nonlinear.relaxation , ...
                 'stagnate'   , false                       , ...
                 'iteration'  , 0                           , ... 
                 'res_history', []                          );

   timer = tic;
   
   converged = false;
   state = state0;

   fprintf('%13s%-26s%-36s\n', '', 'CNV (oil, water)', 'MB (oil, water)');
   

   equation = @(state) equation(tstep, meta.iteration, components, rock, state0, state, dt, bc,dz,p_grad,div,faceConcentrations, system);
   flash=@(fluid) GI_flash(fluid,thermo,options);
   totalFluid=state.totalFluid;
   R=8.3145;
      %%
   % We start with the Newton iterations
   
   while ~ (converged || meta.stopped),
      % Save iteration number in meta info
      meta.iteration = meta.iteration + 1;

      %% 
      % Initial saturation solve
      % 
      % At each Newton step, we start by solving the flash equations and update the liquid
      % saturation variable.
      % %GET RID OF C
      %PROBABLY WONT NEED THIS[C, p] = deal(state.C, state.pressure);
      
%IF tstep=1 and Newton iter =1 then I wanna skip this.So..
if tstep==1 & meta.iteration==1
    state=state0
    p=state.p;
    F=state.F;
    Zi=state.Zi;
    Sw=state.Sw
[p, F, Zi{:}, Sw]=initVariablesADI(state.p, state.F, state.Zi{:}, state.Sw);
%might need to redof lash calc here for ADI purposes
else

    Zi=state.Zi;
    F=state.F;
    p=state.p;
    Sw=state.Sw;
    

    for ji=1:rock.G.cells.num
        totalFluid{ji}.pressure=state.p(ji);
        totalFluid{ji}.mole_fraction=Zi;

        %totalFluid{ji}.mole_fraction=mole_fracs(ji,:);
    end
    
    Xig=[];  state.Xig=[];  %4 components. units=MOLig/MOLg
    Xio=[];  state.Xio=[];%units=MOLio/MOLo
    V=[]; state.V=[];
    Eo=[];  state.Eo=[];%ALREADY Added to each cell
    Eg=[]; state.Eg=[];
    Ew=[];  state.Ew=[];
    Sg=[]; state.Sg=[];
    So=[]; state.So=[];
    rhoLi=[]; state.rhoLi=[];
    rhoGi=[]; state.rhoGi=[];
    rhoL=[];  state.rhoL=[];
    rhoG=[];  state.rhoG=[];

    for j=1:rock.G.cells.num;
    [success_flag,stability_flag,Xiv,Xil,Zgas_vap, Zgas_liq, vapor_frac,cubic_time]=flash(totalFluid{j});
    
    totalFluid{j}.Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
    Xig=[Xig;totalFluid{j}.Xig];
    state.Xig=Xig;
    state.Xig=num2cell(state.Xig,1);

    totalFluid{j}.Xio=Xil(1:3); %units=MOLio/MOLo
    Xio=[Xio;totalFluid{j}.Xio];
    state.Xio=Xio;
    state.Xio=num2cell(state.Xio,1);

    totalFluid{j}.Xwv=Xiv(4); %units=MOLwv/MOLw
    Xwv=[Xwv;totalFluid{j}.Xwv];
    state.Xwv=Xwv;

    totalFluid{j}.Xwl=Xil(4);
    Xwl=[Xwl;totalFluid{j}.Xwl];
    state.Xwl=Xwl;

    totalFluid{j}.V=vapor_frac;
    V=[V;totalFluid{j}.V];
    state.V=V;
    
    totalFluid{j}.Eo=totalFluid{j}.pressure/(Zgas_liq*R*totalFluid{j}.temperature); %ALREADY Added to each cell
    Eo=[Eo;totalFluid{j}.Eo];
    state.Eo=Eo;
    
    totalFluid{j}.Eg=totalFluid{j}.pressure/(Zgas_vap*R*totalFluid{j}.temperature); %ALREADY Added to each cell
    Eg=[Eg;totalFluid{j}.Eg];
    state.Eg=Eg;
    
    totalFluid{j}.Ew=55.5; 
    Ew=[Ew;totalFluid{j}.Ew];
    state.Ew=Ew;
    
    totalFluid{j}.rhoLi=totalFluid{j}.pressure*totalFluid{j}.components.MW(:)/(Zgas_liq*R*totalFluid{j}.temperature);
    rhoLi=[rhoLi;totalFluid{j}.rhoLi];
    state.rhoLi=rhoLi;

    totalFluid{j}.rhoL=sum(totalFluid{j}.rhoLi.*totalFluid{j}.Zi);
    rhoL=[rhoL;totalFluid{j}.rhoL];
    state.rhoL=rhoL;

    totalFluid{j}.rhoGi=totalFluid{j}.pressure*totalFluid{j}.components.MW(i)/(Zgas_vap*R*totalFluid{j}.temperature);
    rhoLi=[rhoGi;totalFluid{j}.rhoGi];
    state.rhoGi=rhoGi;

    totalFluid{j}.rhoG=sum(totalFluid{j}.rhoGi.*totalFluid{j}.Zi);
    rhoG=[rhoG;totalFluid{j}.rhoG];
    state.rhoG=rhoG;
    
    totalFluid{j}.Sg=state.V(j)*state.F(j)/state.Eg(j);
    Sg=[Sg;totalFluid{j}.Sg];
    state.Sg=Sg;
    
    totalFluid{j}.So=1-totalFluid{j}.Sg-state.Sw(j);
    So=[So;totalFluid{j}.So];
    state.So=So;
    state.So(j)=totalFluid{j}.So;
    
   
    end
    state.totalFluid=totalFluid;
    %{
    totalFluid{i}.Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
    totalFluid{i}.Xio=Xil(1:3); %units=MOLio/MOLo
    totalFluid{i}.Xwv=Xiv(4); %units=MOLwv/MOLw
    totalFluid{i}.Xwl=Xil(4);
    totalFluid{i}.V=vapor_frac;
    totalFluid{i}.Eo=totalFluid{i}.pressure/(Zgas_liq*R*totalFluid{i}.fluid.temperature); %ALREADY Added to each cell
    totalFluid{i}.Eg=totalFluid{i}.pressure/(Zgas_vap*R*totalFluid{i}.fluid.temperature); %ALREADY Added to each cell
    totalFluid{i}.Ew=55.5; 
%THIS IS EVERYTHING NOT BEING CHANGED BY NEWTON SOLVER
%}
end

      %%
      % The residual equations for the whole system (pressure, total concentrations,
      % liquid saturation) are assembled.
      %
      
      eqs = equation(state);

      %%
      % We call a standard linear solve to compute the Newton step |dx|.
      
      dx = SolveEqsADI(eqs, []);
      
      %UPDATE STATE WILL CHANGE
      
            %%
      % We update |state|. see below for equation
      %
      state      = updateBWstate(state, dx, 3);
      %STOPPED HERE
      %%
      % We compute the residual values by calling |getResiduals|.
      % This function detects oscillation and stagnation
      
      [meta, residuals] = getResiduals(meta, eqs, system, false);

      %%
      % We test for convergence. Here we use a very stringent test given by a max norm
      %
      fprintf('Newton iteration %d, max residual: %5g\n', meta.iteration, ...
              max(abs(residuals)));
      converged = all(residuals <= system.nonlinear.tol); 

      if meta.stagnate,
         warning('newt:stagnate', 'Non-linear solver stagnated...')
      end
      
      if meta.iteration > system.nonlinear.maxIterations
         meta.stopped = true;
      end
      
   end


   if meta.stopped,
      warning('newt:maxit', ...
              ['Non-linear solver did not converge, stopped ', ...
               'by max iterations...']);
      convergence = false;
   else
      convergence = true;
   end
   
   dispif(mrstVerbose, 'Completed %d iterations in %1.2f s\n', ...
          meta.iteration, toc(timer));

end



      


    
    
    
    
    
    
    
    
    
    
    
    
    
    
