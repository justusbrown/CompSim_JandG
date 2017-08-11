%%function setupBWcontrols

function [bc]=setupBWcontrols(system)
rock=system.rock;
outfluxFluid=system.outfluxFluid;
influxFluid=system.influxFluid;
influx_rate=system.influx_rate;
thermo=system.thermo;
options=system.options;
nComp=system.nComp;
R = thermo.R;

%%
%Setup the dirichlet condition (outflux)
%Note that some of these could be user inputs and are changeable
bound_cell_out=1;
bound_faces = rock.G.cells.faces(rock.G.cells.facePos(bound_cell_out):rock.G.cells.facePos(bound_cell_out + 1) - ...
                              1 , :);
bc.dirichlet.faces =  bound_faces(bound_faces(:, 2) == 2 , 1);

bc.dirichlet.pressure = outfluxFluid.pressure;
bc.dirichlet.fluid=outfluxFluid
%%
% flash for surface properties
[success_flag,stability_flag,Xiv,Xil,Zgas_vap, Zgas_liq, vapor_frac,cubic_time]=GI_flash(bc.dirichlet.fluid,thermo,options);

bc.dirichlet.Xig=Xiv;
bc.dirichlet.Xio=Xil;
bc.dirichlet.So=.25;  %USER INPUT
bc.dirichlet.Sg=.30;  %USER INPUT
bc.dirichlet.Sw=1-bc.dirichlet.So-bc.dirichlet.Sg; %USER INPUT
bc.dirichlet.V=vapor_frac;
bc.dirichlet.Zi=bc.dirichlet.Xig.*bc.dirichlet.V+bc.dirichlet.Xio.*(1-bc.dirichlet.V);
bc.dirichlet.Eo=bc.dirichlet.pressure/(Zgas_liq*R*bc.dirichlet.fluid.Temp); 
bc.dirichlet.Eg=bc.dirichlet.pressure/(Zgas_vap*R*bc.dirichlet.fluid.Temp); 
bc.dirichlet.F=bc.dirichlet.Eo*bc.dirichlet.So+bc.dirichlet.Eg*bc.dirichlet.Sg;
bc.dirichlet.Ew=55.5; %CORRELLATION

%%
%Setup the Neumann condition (influx)
bound_cell_in=30;
bc.in.influx_cells=bound_cell_in;
bc.in.fluid=influxFluid;
bc.in.pressure=influxFluid.pressure;
[success_flag,stability_flag,bc.in.Xiv,bc.in.Xil, bc.in.Zgas_vap, bc.in.Zgas_liq,bc.in.vapor_frac,cubic_time]=GI_flash(bc.in.fluid,thermo,options)

bc.in.Zi=bc.in.Xiv.*bc.in.vapor_frac+bc.in.Xil.*(1-bc.in.vapor_frac);
bc.in.Zi=[bc.in.Zi];

bc.in.Eg=bc.in.pressure/(bc.in.Zgas_vap*R*bc.in.fluid.Temp); 
bc.in.Eo=bc.in.pressure/(bc.in.Zgas_liq*R*bc.in.fluid.Temp); 

for ic=1:nComp
    bc.in.C_influx(ic)=influx_rate*(bc.in.Xiv(ic)*bc.in.vapor_frac+bc.in.Xil(ic)*(1-bc.in.vapor_frac)); %Mols/sec
end

bc.in.T_influx=influx_rate*(bc.in.Eg*bc.in.vapor_frac+bc.in.Eo*(1-bc.in.vapor_frac)); %Mols/sec

bc.in.water_influx=influx_rate/system.mv_water; %Mols/sec: (m^3/s * mol/m^3)
end
