
function [state0]=initTotalFluid_JandG(rock,component,Temp,pressure, options,thermo);

numCells=rock.G.cells.num;
R = getR_JandG();
totalFluid=cell(numCells,1);

for i=1:numCells
    totalFluid{i}=addMixture(component, Temp, pressure);
    
    
[success_flag,stability_flag,Xiv,Xil,Zgas_vap,Zgas_liq,vapor_frac,cubic_time]=GI_flash(totalFluid{i},thermo,options);

totalFluid{i}.Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
totalFluid{i}.Xio=Xil(1:3); %units=MOLio/MOLo
totalFluid{i}.Xwv=Xiv(4); %units=MOLwv/MOLw
totalFluid{i}.Xwl=Xil(4);
totalFluid{i}.V=vapor_frac;
totalFluid{i}.So=.25;
totalFluid{i}.Sg=.30;
totalFluid{i}.Sw=1-totalFluid{i}.So-totalFluid{i}.Sg;
totalFluid{i}.Zi=totalFluid{i}.Xig.*totalFluid{i}.V+totalFluid{i}.Xio.*(1-totalFluid{i}.V); %ALREADY Added to each cell
%totalFluid{i}.Zi=num2cell(totalFluid{i}.Zi,1);
totalFluid{i}.Eo=totalFluid{i}.pressure/(Zgas_liq*R*totalFluid{i}.temperature); %ALREADY Added to each cell
totalFluid{i}.Eg=totalFluid{i}.pressure/(Zgas_vap*R*totalFluid{i}.temperature); %ALREADY Added to each cell
totalFluid{i}.F=(totalFluid{i}.Eo.*totalFluid{i}.So+totalFluid{i}.Eg.*totalFluid{i}.Sg);%ALREADY Added to each cell
totalFluid{i}.Ew=55.5; 
%COULD MAYBE DEFINE STUFF BELOW OUTSIDE LOOP
totalFluid{i}.muL=1e-3;
totalFluid{i}.muG=1e-5;
totalFluid{i}.cl    = 4.4e-5/atm;  % Compressibility
totalFluid{i}.p_ref = 1*atm;       % Reference pressure


end

state.totalFluid=totalFluid;
%{
%ONLY HIGHLIGHTING THE PRIMARIES
   %state.Zi=totalFluid.Zi;
   state.vectorFluid=cell2mat(state.totalFluid);
   state.p=zeros(rock.G.cells.num,1);
   state.F=zeros(rock.G.cells.num,1);
   state.Sw=zeros(rock.G.cells.num,1);
   for iT=1:rock.G.cells.num;
       state.p(iT)=state.totalFluid{iT}.pressure;
       state.F(iT)=state.totalFluid{iT}.F;
       state.Sw(iT)=state.totalFluid{iT}.Sw;
   end %THIS IS ALL TO MAKE INPUTTING THIS INTO SOLVEFI EASY AND EQASSEMBLER
%}
   
state0=state;

   [state0.Xig, state0.Xio, state0.Xwv, state0.Xwl, state0.Eo, state0.Eg, state0.Ew, state0.So, state0.V]=variableCall_JandG(rock,state0);
   state0.Xig=state0.Xig', state0.Xio=state0.Xio', state0.Xwv=state0.Xwv', state0.Xwl=state0.Xwl', state0.Eo=state0.Eo', state0.Eg=state0.Eg', state0.Ew=state0.Ew', state0.So=state0.So', state0.V=state0.V';
    
   [state0.p,state0.F,state0.Sw,state0.Zi]=primeVars_JandG(rock, state0);
   state0.p=state0.p', state0.F=state0.F', state0.Sw=state0.Sw', state0.Zi=state0.Zi';


end
