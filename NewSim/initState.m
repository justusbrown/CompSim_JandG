%NOW INITIALIZING THE STATE
%My initial thoughts on this
%1)Are these just intial values or are updates occuring in these vairables
%2)Some values (So, Sg, and Sw are equal to direchlet?, I'm confused on the
%difference
function [state0]=initState(rock,fluid,options,thermo,Zgas_vap,Zgas_liq);

[success_flag,stability_flag,Xiv,Xil,V,cubic_time]=GI_flash(fluid,thermo,options);
R=getR()
numCells=rock.G.cells.num;

state.Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
state.Xio=Xil(1:3); %units=MOLio/MOLo
state.Xwv=Xiv(4); %units=MOLwv/MOLw
state.Xwl=Xil(4);
state.vapor_frac=vapor_frac;
%CHANGED BELOW 7/19 JB
state.So=.25;
state.Sg=.30;
state.Sw=1-state.So-state.Sg;
%Gage, let me know where you got this formula, I may just be confused
state.Zi=state.Xig*state.Sg+state.Xio*state.So;
%I need change pressure and man up to make sure this is right
state.Eo=state.pressure/(Zgas_liq*R*fluid.temperature); 
state.Eg=state.pressure/(Zgas_vap*R*fluid.temperature); 
state.F=state.Eo*state.So+state.Eg*state.Sg;





state.C=C_initial(1:4).*ones(numCells,1);
state.Xig=vapor_y(1:3).*ones(numCells,1);
state.Xio=liquid_x(1:3).*ones(numCells,1);
state.So=(1-vapor_frac).*ones(numCells,1);
state.Zi=(state.Xig.*vapor_frac+state.Xio.*state.So).*ones(numCells,1);
state.Eo=sum(state.Xio)./rock.pv;
state.Eg=sum(state.Xig)./rock.pv;%HERE I AM TREATING MOLES AND MOLE FRACTIONS AS EQUIVALENT BUT I KNOW THATS NOT OKAY
state.Ew=(liquid_x(4)+vapor_y(4)).*ones(numCells,1)./rock.pv;
state.F=state.Eg.*vapor_frac+state.Eo.*state.So;
state.Sw=vapor_y(4)/sum(vapor_y)+liquid_x(4)/sum(liquid_x);
state.pressure=outflux_p*ones(numCells,1);
state.cwL=liquid_x(4).*ones(numCells,1);
state.cwV=vapor_y(4).*ones(numCells,1);
state.Cw=(vapor_y(4)*vapor_frac+liquid_x(4)*(1-vapor_frac)).*ones(numCells,1); %IM NOT SURE THIS IS ASSEMBLED CORRECTLY, BUT IM TRYING NOT TO FIXATE ON INDIVIDUAL LINES
state.Sg=vapor_frac;%NEED TO KNOW IF VAPOR_FRAC INCLUDES WATER. I THINK IT DOES AND SO SOME THINGS NEED FIXING

state.fluid=initialFluid
%CONSTANT VISCOSITIES & COMPRESSIBILITY (which I don't know if I need atm)
state.fluid.muL=1e-3;
state.fluid.muG=1e-5;
state.fluid.cl    = 4.4e-5/atm  % Compressibility
state.fluid.p_ref = 1*atm       % Reference pressure
