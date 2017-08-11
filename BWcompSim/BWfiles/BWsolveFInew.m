function [state, conv] = BWsolveFInew(tstep, system, ops, state0, bc, ...
          equation);
      
      [components, dt, dz, p_grad, div, faceConcentrations]=deal(system.components,...
         system.options.dt, ops.dz, ops.p_grad, ops.div, ops.faceConcentrations);
   
     [thermo, rock]=deal(system.thermo, system.rock);
     
   meta = struct('converged'  , false                       , ...
                 'stopped'    , false                       , ...
                 'relax'      , system.options.nonlinear.relaxation , ...
                 'stagnate'   , false                       , ...
                 'iteration'  , 0                           , ... 
                 'res_history', []                          );

   timer = tic;
   
   converged = false;
   state = state0;

   fprintf('%13s%-26s%-36s\n', '', 'CNV (oil, water)', 'MB (oil, water)');
   

   equation = @(state) equation(tstep, meta.iteration, state0, state, bc, system, ops);
   flash=@(fluid) GI_flash(fluid,thermo,options);
   totalFluid=state.totalFluid;
   R=thermo.R;
   options=system.options;
   Temp=system.Temp;
   Cells=system.Cells;
   nCell=system.nCell;
   nComp=system.nComp;
   MW=vertcat(components.MW);
   muL=1e-3;
   muG=1e-5;
   cl    = system.cl;  % Compressibility
   p_ref = system.p_ref;
   
   totalFluid=state.totalFluid;
   





      %%
   % We start with the Newton iterations
   
   while ~ (converged || meta.stopped),
      % Save iteration number in meta info
      meta.iteration = meta.iteration + 1;
      
      init_F=state.F;
      init_Zi=state.Zi;
      for i=1:nCell
          totalFluid(i).pressure=state.p(i);
          totalFluid(i).Sw=state.Sw(i);
          totalFluid(i).F=state.F(i);
          totalFluid(i).Zi=[state.Zi{1}(i), state.Zi{2}(i), ...
              state.Zi{3}(i), state.Zi{4}(i), state.Zi{5}(i), state.Zi{6}(i)];
          totalFluid(i).Zi=mynormalize(abs(totalFluid(i).Zi));
      end

      state.totalFluid=totalFluid;
      
      Xig=[]; %units=MOLig/MOLg
      Xio=[]; %units=MOLio/MOLo
      V=[];
      Eo=[];
      Eg=[];
      Ew=[]; 
      rhoL=[];
      rhoG=[];
      Zi=[];
      
      for i=1:nCell
      
      [success_flag,stability_flag,Xiv,Xil, vapor_frac, Zgas_vap,Zgas_liq,cubic_time]=GI_flash(totalFluid(i), thermo,options);
      
        totalFluid(i).Xig=Xiv; %units=MOLig/MOLg
        Xig=[Xig;totalFluid(i).Xig];
        state.Xig=Xig;
        state.Xig=num2cell(state.Xig,1);

        totalFluid(i).Xio=Xil; %units=MOLio/MOLo
        Xio=[Xio;totalFluid(i).Xio];
        state.Xio=Xio;
        state.Xio=num2cell(state.Xio,1);

        totalFluid(i).V=vapor_frac;
        V=[V;totalFluid(i).V];
        state.V=V;
        
        totalFluid(i).Eo=totalFluid(i).pressure/(Zgas_liq*R*totalFluid(i).Temp); %ALREADY Added to each cell
        if isnan(totalFluid(i).Eo)
            totalFluid(i).Eo=0;
        end
        Eo=[Eo;totalFluid(i).Eo];
        state.Eo=Eo;

        totalFluid(i).Eg=totalFluid(i).pressure/(Zgas_vap*R*totalFluid(i).Temp); %ALREADY Added to each cell
        if isnan(totalFluid(i).Eg)
            totalFluid(i).Eg=0;
        end
        Eg=[Eg;totalFluid(i).Eg];
        state.Eg=Eg;
        
        totalFluid(i).rhoL=totalFluid(i).Eo*sum(MW'.*totalFluid(i).Xio);
        rhoL=[rhoL; totalFluid(i).rhoL];
        state.rhoL=rhoL;

        totalFluid(i).rhoG=totalFluid(i).Eg*sum(MW'.*totalFluid(i).Xig);
        rhoG=[rhoG; totalFluid(i).rhoG];
        state.rhoG=rhoG;
        
        totalFluid(i).Zi=totalFluid(i).Xig.*totalFluid(i).V+totalFluid(i).Xio.*(1-totalFluid(i).V); %ALREADY Added to each cell
        Zi=[Zi;totalFluid(i).Zi];
        state.Zi=Zi;
        state.Zi=num2cell(state.Zi,1);
        totalFluid(i).mole_fraction=Zi;
      end
      


      state.totalFluid=totalFluid;
      
      %}
      
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
      state      = updateBWstate(state, dx, system.nComp);
      %STOPPED HERE
      %%
      % We compute the residual values by calling |getResiduals|.
      % This function detects oscillation and stagnation
      
      [meta, residuals] = BWgetResiduals(meta, eqs, system, false);

      %%
      % We test for convergence. Here we use a very stringent test given by a max norm
      %
      fprintf('Newton iteration %d, max residual: %5g\n', meta.iteration, ...
              max(abs(residuals)));
      converged = all(residuals <= system.options.nonlinear.tol); 

      if meta.stagnate,
         warning('newt:stagnate', 'Non-linear solver stagnated...')
      end
      
      if meta.iteration > system.options.nonlinear.maxIterations
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

      
      
      