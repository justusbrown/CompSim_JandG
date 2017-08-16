%% 
% Update function
%
% Given a Newton step, the state variables are updated. We cap the saturation variable.
%
function state = updateBWstate(state, dx, system);

   dp = dx{1};
   dm_w=dx{2};
   for ic = 1 : system.nComp
      dm_i{ic} = dx{2+ic};
   end
   
    
   step = 1;

   state.p = state.p + step*dp;
   state.m_w=state.m_w + step*dm_w;
   for ic = 1 : system.nComp
      state.m_i{ic} = state.m_i{ic} + step*dm_i{ic}; % max(0, state.Zi{ic} + step*dZi{ic});
      %NOT SURE WE NEED MAX HERE
   end
   
   %TO GET RID OF ANY OUTPUTS WHICH ARE NONSENSICAL
for j=1:system.nCell
       
if any(state.p(j)<0)
     state.p(j)=0;
end
   
if any(state.m_w(j)<0)
     state.m_w(j)=0;
end

    for k=1:system.nComp
      
if any(state.m_i{k}(j)<0)
     state.m_i{k}(j)=0;
end

    end
    
end
  
end