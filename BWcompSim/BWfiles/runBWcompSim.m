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
[ops, T]=setupBWsystem(system,bc);
system.rock.T=T;
%%
%Initialize the state
[state0]=initBWstate(system);

%%
%Begin solving the system
%[dt,total_time,steps,t,nComp]=deal(system.dt,system.total_time,...
    %system.steps, system.t, system.nComp);


for tstep = 1 : numel(system.options.steps)
    
      [state, conv] = BWsolveFInew(tstep, system, ops, state0, bc, ...
          @BWeqAssembler);
      
      
       dt = steps(tstep); 
       system.options.dt=dt
       
      if ~(conv)
         error('Convergence failed. Try smaller time steps.')
         return
      end

   if param.do_save
      save(fullfile(param.output_dir, sprintf('state%05d.mat', tstep)), 'state') 
   end
   state0 = state 
   
end


    
    