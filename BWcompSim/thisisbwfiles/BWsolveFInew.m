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
   

   equation = @(state) equation(state0, state, bc, system, ops);
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
      
      init_m_i=state.m_i;
      init_m_w=state.m_w;
      
      for i=1:nCell
          totalFluid(i).pressure=state.p(i);
          totalFluid(i).m_w=state.m_w(i);
          %CAN I MAKE THIS ANOTHER FOR LOOP, NOT THINKING SUPER WELL (I THINK THAT THIS WORKS)?
          %for j=1:nComp
          %totalFluid(i).m_i = state.m_i{j}(i);
          %end
          totalFluid(i).m_i=[state.m_i{1}(i), state.m_i{2}(i), ...
              state.m_i{3}(i), state.m_i{4}(i), state.m_i{5}(i), state.m_i{6}(i)];
      end
      
     
      totalFluid(i).Zi=totalFluid(i).m_i./sum(totalFluid(i).m_i);
      
      state.totalFluid=totalFluid;
      
      Xig=[]; %units=MOLig/MOLg
      Xio=[]; %units=MOLio/MOLo
      V=[];
      Eo=[];
      Eg=[];
      Ew=[]; 
      rhoL=[];
      rhoG=[];
      m_w=[];
      
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

        totalFluid(i).Ew=55.5/system.litre; 
        Ew=[Ew;totalFluid(i).Ew];
        state.Ew=Ew;
        
        totalFluid(i).rhoL=totalFluid(i).Eo*sum(MW'.*totalFluid(i).Xio);
        rhoL=[rhoL; totalFluid(i).rhoL];
        state.rhoL=rhoL;

        totalFluid(i).rhoG=totalFluid(i).Eg*sum(MW'.*totalFluid(i).Xig);
        rhoG=[rhoG; totalFluid(i).rhoG];
        state.rhoG=rhoG;
        
        totalFluid(i).m_w=totalFluid(i).Ew*(1 - sum(totalFluid(i).m_i)*(1-totalFluid(i).V)/totalFluid(i).Eo - sum(totalFluid(i).m_i)*totalFluid(i).V/totalFluid(i).Eg);
        m_w=[m_w; totalFluid(i).m_w];
        state.m_w=m_w;
        
       end
       
       state.totalFluid=totalFluid;  
       eqs = equation(state);    
       dx = SolveEqsADI(eqs, []);
       state      = updateBWstate(state, dx, system.nComp);
       
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
       
       




      
      
          
          
          
          
      
     
      
         
      