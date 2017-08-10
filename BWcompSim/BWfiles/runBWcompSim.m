%%
mrstModule add ad-fi 
mrstModule add ad-core

gravity off

%%
%See function inputData() to change user inputs
[system]=BWinputData()

%%
%Setup the controls/bounds
bc=setupBWcontrols(system);

%%
%Setup the system/operators
ops=setupBWsystem(system,bc);

%%
%Initialize the state
[state0]=initBWstate(system);

%%
%Begin solving the system
[dt,total_time,steps,t,nComp]=deal(system.options.dt,system.options.total_time,...
    system.options.steps, system.options.t, system.nComp);


for tstep = 1 : numel(steps)
    
      [state, conv] = BWsolveFI(tstep, system, ops, thermo, rock, state0, bc, ...
          @BWeqAssembler, options);
      
      
       dt = steps(tstep); 
       
      if ~(conv)
         error('Convergence failed. Try smaller time steps.')
         return
      end

   if param.do_save
      save(fullfile(param.output_dir, sprintf('state%05d.mat', tstep)), 'state') 
   end
   state0 = state 
   
end


    
    