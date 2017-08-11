function [state0]=initBWstate(system);

thermo=system.thermo;
options=system.options;
rock=system.rock;
components=system.components;
MW=vertcat(components.MW);
Temp=system.Temp;
Cells=system.Cells;
nCell=system.nCell;
nComp=system.nComp;
R=thermo.R;
totalFluid=system.totalFluid;
%pressure=totalFluid.pressure;
muL=1e-3;
muG=1e-5;
cl    = system.cl;  % Compressibility
p_ref = system.p_ref;       % Reference pressure

%for i=1:nCell
 %   p(i)=totalFluid(i).pressure;
  %  Zi(i,:)=totalFluid(i).mole_fraction;
%end
%p=p';
%Zi=num2cell(Zi, 1);

Zi=[];
Xig=[]; %4 componentss. units=MOLig/MOLg
Xio=[]; %units=MOLio/MOLo
V=[];
So=[];
Sg=[];
Sw=[];
Eo=[]; 
Eg=[]; 
F=[];
Ew=[]; 
%rhoLi=[];
%rhoGi=[];
rhoL=[];
rhoG=[];

%pressure=cat(1,totalFluid.pressure);

for i=1:nCell
totalFluid(i).Zi=totalFluid(i).mole_fraction;
end

for i=1:nCell    
    
[success_flag,stability_flag,Xiv,Xil,Zgas_vap,Zgas_liq,vapor_frac,cubic_time]=GI_flash(totalFluid(i), thermo,options);

totalFluid(i).Xig=Xiv; %4 components. units=MOLig/MOLg
Xig=[Xig;totalFluid(i).Xig];
state.Xig=Xig;
state.Xig=num2cell(state.Xig,1);

totalFluid(i).Xio=Xil; %units=MOLio/MOLo
Xio=[Xio;totalFluid(i).Xio];
state.Xio=Xio;
state.Xio=num2cell(state.Xio,1);

totalFluid(i).V=vapor_frac;
V=[V;totalFluid(i).V];
state.V=V;

totalFluid(i).So=.25;
So=[So;totalFluid(i).So];
state.So=So;

totalFluid(i).Sg=.30;
Sg=[Sg;totalFluid(i).Sg];
state.Sg=Sg;

totalFluid(i).Sw=1-totalFluid(i).So-totalFluid(i).Sg;
Sw=[Sw;totalFluid(i).Sw];
state.Sw=Sw;

totalFluid(i).Zi=totalFluid(i).Xig.*totalFluid(i).V+totalFluid(i).Xio.*(1-totalFluid(i).V); %ALREADY Added to each cell
Zi=[Zi;totalFluid(i).Zi];
state.Zi=Zi;
state.Zi=num2cell(state.Zi,1);
totalFluid(i).mole_fraction=Zi;

totalFluid(i).Eo=totalFluid(i).pressure/(Zgas_liq*R*totalFluid(i).Temp); %ALREADY Added to each cell
Eo=[Eo;totalFluid(i).Eo];
state.Eo=Eo;

totalFluid(i).Eg=totalFluid(i).pressure/(Zgas_vap*R*totalFluid(i).Temp); %ALREADY Added to each cell
Eg=[Eg;totalFluid(i).Eg];
state.Eg=Eg;

totalFluid(i).rhoL=totalFluid(i).Eo*sum(MW'.*totalFluid(i).Xio);
rhoL=[rhoL; totalFluid(i).rhoL];
state.rhoL=rhoL;

totalFluid(i).rhoG=totalFluid(i).Eg*sum(MW'.*totalFluid(i).Xig);
rhoG=[rhoG; totalFluid(i).rhoG];
state.rhoG=rhoG;

%totalFluid(i).rhoLi=totalFluid(i).pressure.*MW/(Zgas_liq*R*totalFluid(i).Temp);
%rhoLi=[rhoLi;totalFluid(i).rhoLi];
%state.rhoLi=rhoLi;

%totalFluid(i).rhoL=totalFluid(i).rhoLi.*totalFluid(i).Zi;
%totalFluid(i).rhoL=sum(totalFluid(i).rhoL);
%rhoL=[rhoL;totalFluid(i).rhoL];
%state.rhoL=rhoL;
%Not summing properly

%totalFluid(i).rhoGi=totalFluid(i).pressure.*MW/(Zgas_vap*R*totalFluid(i).Temp);
%rhoLi=[rhoGi;totalFluid(i).rhoGi];
%state.rhoGi=rhoGi;

%totalFluid(i).rhoG=totalFluid(i).rhoGi.*totalFluid(i).Zi;
%totalFluid(i).rhoG=sum(totalFluid(i).rhoG);
%rhoG=[rhoG;totalFluid(i).rhoG];
%state.rhoG=rhoG;
%Not summing properly

totalFluid(i).F=(totalFluid(i).Eo.*totalFluid(i).So+totalFluid(i).Eg.*totalFluid(i).Sg);%ALREADY Added to each cell
F=[F;totalFluid(i).F];
state.F=F;

totalFluid(i).Ew=55.5; 
Ew=[Ew;totalFluid(i).Ew];
state.Ew=Ew;


end



state.totalFluid=totalFluid;

for i=1:nCell
   state.p(i)=totalFluid(i).pressure;
end
state.p=state.p'

  % state.p=zeros(rock.G.cells.num,1);
   
  % for iT=1:rock.G.cells.num;
       %state.p(iT)=state.totalFluid{iT}.pressure;
   %end
   %state.p=state.p'
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

   %[state0.Xig, state0.Xio, state0.Xwv, state0.Xwl, state0.Eo, state0.Eg, state0.Ew, state0.So, state0.V]=variableCall_JandG(rock,state0);
   %state0.Xig=state0.Xig', state0.Xio=state0.Xio', state0.Xwv=state0.Xwv', state0.Xwl=state0.Xwl', state0.Eo=state0.Eo', state0.Eg=state0.Eg', state0.Ew=state0.Ew', state0.So=state0.So', state0.V=state0.V';
    
   %[state0.p,state0.F,state0.Sw,state0.Zi]=primeVars_JandG(rock, state0);
   %state0.p=state0.p', state0.F=state0.F', state0.Sw=state0.Sw', state0.Zi=state0.Zi';


end
