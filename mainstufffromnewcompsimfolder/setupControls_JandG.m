%%function setupControls_JandG

function [bc]=setupBWcontrols(rock,outfluxFluid,influxFluid,influx_rate, thermo,options, system)

R = getR_JandG();

%%
%Setup the dirichlet condition
%Note that some of these could be user inputs and are changeable
bound_cell_out=1;
bound_faces = rock.G.cells.faces(rock.G.cells.facePos(bound_cell_out):rock.G.cells.facePos(bound_cell_out + 1) - ...
                              1 , :);
bc.dirichlet.faces =  bound_faces(bound_faces(:, 2) == 2 , 1);

bc.dirichlet.pressure = outfluxFluid.pressure;
bc.dirichlet.fluid=outfluxFluid
% flash for surface properties
[success_flag,stability_flag,Xiv,Xil,Zgas_vap, Zgas_liq, vapor_frac,cubic_time]=GI_flash(bc.dirichlet.fluid,thermo,options)
bc.dirichlet.Xig=Xiv;
bc.dirichlet.Xio=Xil;
%I AM TAKING WATER OUT OF FLASH
%bc.dirichlet.Xwv=Xiv;
%bc.dirichlet.Xwl=Xil(4);
%SINCE THIS IS A BOUNDARY CONDITION, I GUESS IT CAN BE A SINGLE VALUE? JB
%7/21
%THESE CAN BE MOVED TO USER INPUTS inputData
bc.dirichlet.So=.25; 
bc.dirichlet.Sg=.30; 
bc.dirichlet.Sw=1-bc.dirichlet.So-bc.dirichlet.Sg;
bc.dirichlet.V=vapor_frac;
bc.dirichlet.Zi=bc.dirichlet.Xig.*bc.dirichlet.V+bc.dirichlet.Xio.*(1-bc.dirichlet.V); %CHANGED V_FRAC DOESN'T EQUAL SG %CHANGED!!! jb-07/18
bc.dirichlet.Eo=bc.dirichlet.pressure/(Zgas_liq*R*bc.dirichlet.fluid.temperature); 
bc.dirichlet.Eg=bc.dirichlet.pressure/(Zgas_vap*R*bc.dirichlet.fluid.temperature); 
bc.dirichlet.F=bc.dirichlet.Eo*bc.dirichlet.So+bc.dirichlet.Eg*bc.dirichlet.Sg;
bc.dirichlet.Ew=55.5;

%{
[bc.dirichlet.Cw,bc.dirichlet.cwg,bc.dirichlet.cwl]=computeWaterComp_JandG(bc.dirichlet.pressure,bc.dirichlet.Zi,bc.dirichlet.Xio,bc.dirichlet.V,system);
%}
%SOMETHING IS MAJORLY WRONG WITH THESE WATER NUMBERS!!!

%bc.dirichlet.Ew=? DONT KNOW IF NEEDED HERE BUT DO NEED FOR INITIAL STATE
%gr-7/19%YES, WE WILL NEED INITIAL VALUES FOR ALL VARIABLE (BOTH PRIMARY AND
%SECONDARY) JB 7/21

%WATERbc.dirichlet.cwL=Xil(4);%SLIGHTLY CONFUSED ON WHAT cwL and Cw is, I know you told me
%WATERbc.dirichlet.cwV=Xiv(4);
%WATERbc.dirichlet.Cw=Xiv(4)*vapor_frac+Xil(4)*(1-vapor_frac) %IM NOT SURE THIS IS ASSEMBLED CORRECTLY, BUT IM TRYING NOT TO FIXATE ON INDIVIDUAL LINES


bound_cell_in=30;
bc.in.influx_cells=bound_cell_in;
[success_flag,stability_flag,bc.in.Xiv,bc.in.Xil,bc.in.vapor_frac,cubic_time]=GI_flash(influxFluid,thermo,options) %DONT SEEM TO USE ANYTHING FROM THIS FLASH TEST...WHY?
for ic=1:3 %DOING WATER SEPERATELY
    bc.in.C_influx(ic)=influx_rate*(bc.in.Xiv(ic)*bc.in.vapor_frac+bc.in.Xil(ic)*(1-bc.in.vapor_frac)); %IM NOT SURE IF THIS IS CORRECT EITHER
end
bc.in.Zi=bc.in.Xiv.*bc.in.vapor_frac+bc.in.Xil.*(1-bc.in.vapor_frac);

[bc.in.Cw,bc.in.cwg,bc.in.cwl]=computeWaterComp_JandG(influxFluid.pressure,bc.in.Zi,bc.in.Xil,bc.in.vapor_frac,system);
bc.in.water_influx=influx_rate*(bc.in.Xil(4)*(1-bc.in.vapor_frac)+bc.in.Xiv(4)*bc.in.vapor_frac);
end
