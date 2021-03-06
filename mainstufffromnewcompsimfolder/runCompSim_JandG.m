 %%
 %MAIN FILE OF CompSim_JandG REPOSITORY
 
 
%IMPORTANT NOTE GR 07/20:
%BravoDome has everything for the solver formated ass Cell arrays using {}.
%We might have to change the character types on these. Shouldn't be
%timely if so.
mrstModule add ad-fi 
mrstModule add ad-core

%%
%See function inputData() to change user inputs
[G,rock,options,thermo,influxFluid,outfluxFluid,initialFluid, influx_rate, system]=inputData()

% Disable gravity
gravity off

%clf
%plotCellData(G,G.cells.indexMap), view(30,50), xlabel('x'), ylabel('y'),zlabel('z'), colorbar

%%
%THIS IS WHERE setupControls WAS gr-07/19
bc=setupControls_JandG(rock,outfluxFluid,influxFluid,influx_rate,thermo,options, system);

%%
%SETUP SYSTEM
   cf = rock.G.cells.faces(:,1);
   nf = rock.G.faces.num;
   nc = rock.G.cells.num;
   %NOTE using rock.pv and rock.poro, don't want an s object like bravodome
   %COMPUTE FULLTRANSMITIBILITY
   rock.Tfull=1 ./ accumarray(cf,1./rock.T,[nf,1]);
   rock.T=rock.Tfull;
   
   %SETUPDISCRETE DIVERGENCE OPERATOR
   [C1,C2,C_div,div]=divOp_JandG(rock.G, cf, nf, nc);
   
   %GAGE 7/19 STARTED SETTING UP SYSTEM. FINSHED DIVERGENCE OPERATOR. 
   %CODE WORKS UP TO THIS POINT. NEED GRADIENT OPERATOR STILL.
   %DECIDED TO KEEP SETUP SYSTEM IN THE CORE FILE. THIS WILL NOT BE THE
   %CASEWITH ASSEMBLING THE EQUATIONS.
   %NEXT STEP IS SETTING UP GRADIENT OPERATOR
   
      % Set up the discrete gradient operator, |grad|. We compute the differences of cell values
   % across each face. It is a linear mapping from cells' to faces' values.
   N = double(rock.G.faces.neighbors);
   index = 1:nf;
   interior = prod(N, 2)~=0;

   C1_interior = sparse(index(interior), N(interior, 1), ones(nnz(interior), 1), nf, nc);
   C2_interior = sparse(index(interior), N(interior, 2), ones(nnz(interior), 1), nf, nc);

   
   %%
   % Compute the boundary contribution to the gradient operator. They corresponds to the
   % external faces where Dirichlet conditions are given. We are careful to use the
   % correct signs.
   % 

   is_dirichlet_faces1 = N(bc.dirichlet.faces, 1) ~= 0;
   is_dirichlet_faces2 = ~is_dirichlet_faces1;

   dirichlet_faces1 = bc.dirichlet.faces(is_dirichlet_faces1);
   dirichlet_faces2 = bc.dirichlet.faces(is_dirichlet_faces2);

   C1_exterior = sparse(index(dirichlet_faces1), ...
                        N(dirichlet_faces1, 1), ...
                        ones(numel(dirichlet_faces1), 1), nf, nc);
   C2_exterior = sparse(index(dirichlet_faces2), ...
                        N(dirichlet_faces2, 2), ...
                        ones(numel(dirichlet_faces2), 1), nf, nc);
                    
   %%
   % The gradient operator is the sum of internal and boundary contributions.
   %
   
   C = C1_interior + C1_exterior - (C2_interior + C2_exterior);

   pressure_bc = sparse(nf, 1);
   pressure_bc(dirichlet_faces1) = - bc.dirichlet.pressure(is_dirichlet_faces1);
   pressure_bc(dirichlet_faces2) = + bc.dirichlet.pressure(is_dirichlet_faces2);

   p_grad = @(p)(C*p + pressure_bc);  

   grad = @(val, bc_val)(grad_JandG(val, bc_val, nf, C, ...
                                is_dirichlet_faces1, dirichlet_faces1, ... 
                                is_dirichlet_faces2, dirichlet_faces2));

%%
   % Set up the gravity term.
   %
   
   z = rock.G.cells.centroids(:, 3);
   fz = rock.G.faces.centroids(:, 3);
   dz = grad(z, fz);
   
      %% 
   % Define function to compute concentrations on faces using cells' values and upwind
   % directions. Here, the _nature_ of the boundary faces (Dirichlet or note) are not
   % allowed to change, although the _values_ on these faces may change.
   %
   
   faceConcentrations = @(flag, conc_c, bc_conc) ...
       faceConcentrations_JandG(flag, conc_c, bc_conc, N, interior, dirichlet_faces2, ...
                          dirichlet_faces1,  bc, nf, nc);
  %THIS IS THE END OF SETUP SYSTEM. I BELIEVE IT IS COMPLETE, and KNOW 
  %THAT EVERYTHING STILL FUNCTIONS.
  
   %SETUP NONLINEAR SOLVER WILL START HERE
   maxIteration=50;
   system.nonlinear=setNonlinearSolverParameters_JandG(maxIteration);
   system.cellwise  = 1:5  % Used in function getResiduals which checks convergence.

   %Done setting up nonlinear solver gr 07/20
   
   %%NOW WE WILL INITIALIZE THE STATE USING initState_JandG
   %gr 07/20
   %STUPID QUESTION, SHOULD WE BE CALLING initialFluid in the
   %initState_JandG sheet as well: [state0]=initState_JandG(rock,fluid,options,thermo);
   [state0]=initTotalFluid_JandG(rock,system.components, system.Temp, outfluxFluid.pressure, options, thermo);
   %MIGHT WANNA GET TOTALFLUID HERE AS WELL 
   %AND THE STATE IS INITIALIZED AS state0



%%
%%LETS START SOLVING THE SYSTEM
dt = 200*day;              % Time step
total_time = 100000*day;   % Total time

steps      = dt*ones(floor(total_time/dt), 1); 
t          = cumsum(steps); 

nComp_C=3; %# OF NON WATER COMPONENTS
%% 
% HERE ARE THE TIME STEPS STARTINg FOR THE ITERATIONS

for tstep = 1 : numel(steps)
   % Call non-linear solver perlAddLink(solvefi)
  [state, conv] = solvefi_JandG_underConstruction(system.components, tstep, system, rock, state0, dt, bc, dz, p_grad, div, faceConcentrations, @eqAssembler_JandG, options);
  %SOLVEFI IS NOT GONNA CONTAIN SYSTEM OR PARAM. IT WILL CONTAIN OTHER
  %THINGS


   dt = steps(tstep); 
   
%EQUATION ASSEMBLING HAS BEEN MOVED TO eqAssembler_JandG.
   %THIS LINE IS JUST A TEST TO SEE IF THERE ARE ANY ERRORS MADE SO FAR: GR
   %07/20
  % TESTLINE:equation=@(state) eqAssembler_JandG(rock,state0,state,dt,bc,dz,p_grad,div,faceConcentrations); 
%EVERYTHING IS WORKING SOMEHOW... maybe not correctly, but its working

   if ~(conv)
      error('Convergence failed. Try smaller time steps.')
      return
   end

   if param.do_save
      save(fullfile(param.output_dir, sprintf('state%05d.mat', tstep)), 'state') 
   end
   state0 = state 
   
end

  





