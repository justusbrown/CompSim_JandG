%%
mrstModule add ad-fi 
mrstModule add ad-core

gravity off

%%
%See function inputData() to change user inputs
[G,rock,options,thermo,influxFluid,outfluxFluid,initialFluid, influx_rate, system]=BWinputData()

%%
%Setup the controls/bounds
bc=setupBWcontrols(rock,outfluxFluid,influxFluid,influx_rate,thermo,options, system);

%%
%Setup the system/operators
ops=setupBWsystem(rock,bc);

%%
%Initialize the state
[state0]=initBWstate(rock,system, initialFluid.pressure, options, thermo);

%%
%Begin solving the system
[dt,total_time,steps,t,nComp]=deal(system.dt,system.total_time,...
    system.steps, system.t, system.nComp);


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


    
    