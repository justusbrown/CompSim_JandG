%% 
% Update function
%
% Given a Newton step, the state variables are updated. We cap the saturation variable.
%
function state = updateState_JandG(state, dx, nComp);

   dp = dx{1};
   dF=dx{2};
   for ic = 1 : nComp
      dZi{ic} = dx{2+ic};
   end
   dSw = dx{nComp + 3};
   
    
   step = 1;

   state.pressure = state.pressure + step*dp;
   state.F=state.F+step*dF;
   for ic = 1 : nComp
      state.Zi{ic} = max(0, state.Zi{ic} + step*dZi{ic});
      %NOT SURE WE NEED MAX HERE
   end
   state.Sw = state.Sw + step*dSw;
   %GOT RID OF min FUNCTION STUFF BUT MIGHT NEED TO GO BACK.
 
   
end
