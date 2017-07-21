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
function [state, convergence] = solvefi_J&G(rock,state0, dt, bc, dz,p_grad,div,faceConcentrations, equation,options,  varargin)
%
   opt = struct('verbose', false);
   opt = merge_options(opt, varargin{:});

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
   
   %DON'T REALLY UNDERSTAND THIS PART
   %LITERALLY JUST MAKING IT SO THAT WE TYPE EQUATION INSTEAD @eqsAssemble
  % ANY "SYSTEM" HAS TO GO gr 07/20
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
[Sw,p]=deal(state.Sw,state.pressure);
      %JUST REALIZED, ITHINK Sg NEEDS TO DEPEND ON Sw AND NOT VICE VERSA
state.fluid.pressure=p;
[success_flag,stability_flag,Xiv,Xil,Zgas_vap, Zgas_liq, vapor_frac,cubic_time]=flash(state.fluid);
state.Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
state.Xio=Xil(1:3); %units=MOLio/MOLo
state.Xwv=Xiv(4); %units=MOLwv/MOLw
state.Xwl=Xil(4);
state.V=vapor_frac;
%CHANGED BELOW 7/19 JB
%state.So=.25; THIS WILL BE SOLVED FOR
%state.Sg=.30; THIS ALSO
%THESE SHOULD BE SOLVED FOR THIS TIME HERE I THINK.
state.Sw=1-state.So-state.Sg;
state.Zi=state.Xig.*state.V+state.Xio.*(1-state.V); 
state.Eo=state.pressure/(Zgas_liq*R*state.fluid.temperature); 
state.Eg=state.pressure/(Zgas_vap*R*state.fluid.temperature); 
state.F=(state.Eo.*state.So+state.Eg.*state.Sg);
%IM TIRED AND CONFUSED LOL, WE ARE GOING TO SOLVE FOR THESE DIFFERENTLY
%THAN I CURRENTLY HAVE IT THOUGH PROBABLY. BUT THIS IS IN THE RIGHT
%DIRECTION GR 07/20
%Still need Ew
      %nc = numel(p);
      %state.sL=state.So+state.Sg; %ADDED!!! jb 7/20
      %init_sL = state.sL; 
      %Cvec = cell2mat(C);
      %[cg, cl, cw, s, Cw] = flash_calculation(Cvec, p, system, init_sL);
      %state.sL = s(:, 2);
      
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
      state      = updateState(state, dx, system);
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

%% 
% Update function
%
% Given a Newton step, the state variables are updated. We cap the saturation variable.
%

function state = updateState(state, dx, system)

   nComp = system.nComp;
   dp = dx{1};
   for ic = 1 : nComp
      dC{ic} = dx{ic + 1};
   end
   dsL = dx{nComp + 2};
    
   step = 1;

   state.pressure = state.pressure + step*dp;
   for ic = 1 : nComp
      state.C{ic} = max(0, state.C{ic} + step*dC{ic});
   end
   state.sL = min(1, max(0, state.sL + step*dsL));
   
end
