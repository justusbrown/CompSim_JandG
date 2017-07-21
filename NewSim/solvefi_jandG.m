%%
%
% perlOneLineDescription(Single step non-linear solver)
% 
% Compute the next |state| of the system for given previous state |state0| and time step
% |dt|.
%
% The discretized residual equation are given by |equation| which is assembled in the
% function |equationCompositional|. 
%%
function [state, convergence] = solvefi_JandG(system, rock,state0, dt, bc, dz,p_grad,div,faceConcentrations, equation,options)
%
   %opt = struct('verbose', false);
   %opt = merge_options(opt, varargin{:});

   fluid = state.fluid;%changed to state, didn't see system anywhere else in our code
   
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
   

   equation = @(state) equation(rock, state0, state, dt, bc,dz,p_grad,div,faceConcentrations); system);
   flash=@(state.fluid) GI_flash(state.fluid,thermo,options);
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
[Sw, p, F, Zi]=deal(state.Sw, state.pressure, state.F, state.Zi);
      %JUST REALIZED, ITHINK Sg NEEDS TO DEPEND ON Sw AND NOT VICE VERSA
state.fluid.pressure=p; %Update the fluids pressure since its seperate from P 
state.fluid.mole_fraction=Zi;%Update fluids mole fraction since its seperate from Zi
[success_flag,stability_flag,Xiv,Xil,Zgas_vap, Zgas_liq, vapor_frac,cubic_time]=flash(state.fluid);

state.Xig=Xiv(1:3);
state.Xio=Xil(1:3); %units=MOLio/MOLo
state.Xwv=Xiv(4); %units=MOLwv/MOLw
state.Xwl=Xil(4);
state.V=vapor_frac;
state.Eo=state.pressure/(Zgas_liq*R*state.fluid.temperature); 
state.Eg=state.pressure/(Zgas_vap*R*state.fluid.temperature); 
state.Ew=55.5; %THIS IS A PLACEHOLDER


%%
%THIS SECTION NEEDS TO BE CHECKED TO SEE IF YOU AGREE
%THE LINE BENEATH THIS I THINK IS NOT NEEDED, SW IS A PRIMARY VARIABLE
%WHICH IS SOLVED FOR IN THE NEWTON ITERATION, WE ONLY NEED THE INITIAL CONDITION TO GET US STARTED JB 7/21
%state.Sw=1-state.So-state.Sg;

%THESE NEED TO BE UPDATED THROUGH PVT INFO JB 7/21


%%
%Still need Ew!!!!!!!!
%MY THOUGHTS: JB 7/21
%MOLARITY IN MY THOUGHTS IS MORE OF A MIXTURE PROPERTY, IF WATER IS A PURE
%SUBSTANCE, DOES THE MOLARITY VARY MUCH? IF IS DOESN'T, WE CAN HARD CODE A
%VALUE (THIS SEEMS TO EASY), IF NOT, THE EW WILL CHANGE AND NEED UPDATING
%THROUGH PVT INFORMATION JB 7/21
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


