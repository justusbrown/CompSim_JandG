function [Xig,Xio,Xwv,Xwl,Eo,Eg,Ew,So,V]=variableCall_JandG(rock,state)

   for iT=1:rock.G.cells.num;
       Xig{iT}=num2cell(state.totalFluid{iT}.Xig, 1);
       Xio{iT}=num2cell(state.totalFluid{iT}.Xio, 1);
       Xwv(iT)=state.totalFluid{iT}.Xwv;
       Xwl(iT)=state.totalFluid{iT}.Xwl;
       Eo(iT)=state.totalFluid{iT}.Eo;
       Eg(iT)=state.totalFluid{iT}.Eg;
       Ew(iT)=state.totalFluid{iT}.Ew;
       So(iT)=state.totalFluid{iT}.So;
       V(iT)=state.totalFluid{iT}.V;
   end
end