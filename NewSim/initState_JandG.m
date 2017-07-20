%NOW INITIALIZING THE STATE
%My initial thoughts on this
%1)Are these just intial values or are updates occuring in these vairables
%2)Some values (So, Sg, and Sw are equal to direchlet?, I'm confused on the
%difference
function [state0]=initState_JandG(rock,fluid,options,thermo);

[success_flag,stability_flag,Xiv,Xil,Zgas_vap,Zgas_liq,vapor_frac,cubic_time]=GI_flash(fluid,thermo,options);
R = getR_JandG();
numCells=rock.G.cells.num;

state.pressure=fluid.pressure;
%I am including state.fluid=initial fluid. this might be left out. same
%for bc
state.fluid=fluid;
state.Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
state.Xio=Xil(1:3); %units=MOLio/MOLo
state.Xwv=Xiv(4); %units=MOLwv/MOLw
state.Xwl=Xil(4);
state.V=vapor_frac;
%CHANGED BELOW 7/19 JB
state.So=.25;
state.Sg=.30;
state.Sw=1-state.So-state.Sg;
%Gage, let me know where you got this formula, I may just be confused
%Thanks Justus, its fixed. gr 07/20
state.Zi=state.Xig*state.V+state.Xio*(1-state.V);
%I need change pressure and man up to make sure this is right
state.Eo=state.pressure/(Zgas_liq*R*state.fluid.temperature); 
state.Eg=state.pressure/(Zgas_vap*R*state.fluid.temperature); 
state.F=state.Eo*state.So+state.Eg*state.Sg;


%CONSTANT VISCOSITIES & COMPRESSIBILITY (which I don't know if I need atm)
state.fluid.muL=1e-3;
state.fluid.muG=1e-5;
state.fluid.cl    = 4.4e-5/atm;  % Compressibility
state.fluid.p_ref = 1*atm;       % Reference pressure

state0=state;

end
