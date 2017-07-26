function [state, convergence] = solvefi_JandG(component, tstep, system, rock,state0, dt, bc, dz,p_grad,div,faceConcentrations, equation,options)
%
   %opt = struct('verbose', false);
   %opt = merge_options(opt, varargin{:});

   
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
   
   thermo=addThermo();
   thermo.EOS=@PREOS;
   equation = @(state) equation(component, rock, state0, state, dt, bc,dz,p_grad,div,faceConcentrations);
   flash=@(fluid) GI_flash(fluid,thermo,options);
   totalFluid=state.totalFluid;
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
else
    %NEED TO REDFINE SOMETHINGS HERE IN UPDATE STAATE
    for j=1:rock.G.cells.num;
    [success_flag,stability_flag,Xiv,Xil,Zgas_vap, Zgas_liq, vapor_frac,cubic_time]=flash(totalFluid{j});
    
    totalFluid{i}.Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
    totalFluid{i}.Xio=Xil(1:3); %units=MOLio/MOLo
    totalFluid{i}.Xwv=Xiv(4); %units=MOLwv/MOLw
    totalFluid{i}.Xwl=Xil(4);
    totalFluid{i}.V=vapor_frac;
    totalFluid{i}.Eo=totalFluid{i}.pressure/(Zgas_liq*R*totalFluid{i}.fluid.temperature); %ALREADY Added to each cell
    totalFluid{i}.Eg=totalFluid{i}.pressure/(Zgas_vap*R*totalFluid{i}.fluid.temperature); %ALREADY Added to each cell
    totalFluid{i}.Ew=55.5; 
%THIS IS EVERYTHING NOT BEING CHANGED BY NEWTON SOLVER

    end
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
      state      = updateState_JandG(state, dx, 3);
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



      


    
    
    
    
    
    
    
    
    
    
    
    
    
    
