%% 
% Update function
%
% Given a Newton step, the state variables are updated. We cap the saturation variable.
%
function state = updateBWstate(state, dx, nComp);

   dp = dx{1};
   dm_w=dx{2};
   for ic = 1 : nComp
      dm_i{ic} = dx{2+ic};
   end
   
    
   step = 1;

   state.p = state.p + step*dp;
   state.m_w=state.m_w + step*dm_w;
   for ic = 1 : nComp
      state.m_i{ic} = state.m_i{ic} + step*dm_i{ic}; % max(0, state.Zi{ic} + step*dZi{ic});
      %NOT SURE WE NEED MAX HERE
   end
 
   
end