%NOW INITIALIZING THE STATE
function [state0]=initState(rock,fluid,options,thermo);

[success_flag,stability_flag,Xiv,Xil,V,cubic_time]=GI_flash(fluid,thermo,options);

numCells=rock.G.cells.num;

state.Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
state.Xio=Xil(1:3); %units=MOLio/MOLo
state.Xwv=Xiv(4); %units=MOLwv/MOLw
state.Xwl=Xil(4);
state.V=V




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
