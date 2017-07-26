
function [p,F,Sw,Zi]=primeVars_JandG(rock, state)
   %p=zeros(rock.G.cells.num,1);
   %F=zeros(rock.G.cells.num,1);
   %Sw=zeros(rock.G.cells.num,1);
   for iT=1:rock.G.cells.num;
       p(iT)=state.totalFluid{iT}.pressure;
       F(iT)=state.totalFluid{iT}.F;
       Sw(iT)=state.totalFluid{iT}.Sw;
       Zi{iT}=num2cell(state.totalFluid{iT}.Zi,1);
   end 
end
