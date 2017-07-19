%%function setupControls_JandG

function [bc]=setupControls_JandG(rock,outfluxFluid,thermo,options)

%CONFUSED ON THIS BLOCK. COULD JUST SAY bc.dirichlet.faces=2
bound_cell_out=1;
bound_faces = G.cells.faces(G.cells.facePos(bound_cell_out):G.cells.facePos(bound_cell_out + 1) - ...
                              1 , :);
bc.dirichlet.faces =  bound_faces(bound_faces(:, 2) == 2 , 1);
%%
%%
%SET UP DIRICHLET CONDITION, CHANGED jb-07/18!!!!
bc.dirichlet.pressure = outfluxFluid.pressure;
bc.dirichlet.fluid=outfluxFluid
% flash for surface properties
[success_flag,stability_flag,Xiv,Xil,vapor_frac,cubic_time]=GI_flash(bc.dirichlet.fluid,thermo,options)
bc.dirichlet.Xig=Xiv(1:3);
bc.dirichlet.Xio=Xil(1:3);
bc.dirichlet.So=.25; %Should be input from user %CHANGED!!! jb-07/18
bc.dirichlet.Sg=.30; %Should be input from user %CHANGED!!! jb-07/18
bc.dirichlet.Sw=1-bc.dirichlet.So-bc.dirichlet.Sg; %for simplicity %CHANGED!!! jb-07/18
bc.dirichlet.Zi=bc.dirichlet.Xig*bc.dirichlet.Sg+bc.dirichlet.Xio*bc.dirichlet.So; %CHANGED V_FRAC DOESN'T EQUAL SG %CHANGED!!! jb-07/18
%THE FOLLOWING NEEDS TO BE FROM PREOS, NEEDS TO BE CHANGED, RATIO OF NUM OF
%MOLS TO GAS VOLUME %CHANGED!!!! jb-07/18
%bc.dirichlet.Eio= 
%bc.dirichlet.Eig=
bc.dirichlet.F=bc.dirichlet.Eo*bc.dirichlet.So+bc.dirichlet.Eg*bc.dirichlet.Sg;

bc.dirichlet.cwL=Xil(4);%SLIGHTLY CONFUSED ON WHAT cwL and Cw is, I know you told me
bc.dirichlet.cwV=Xiv(4);
bc.dirichlet.Cw=Xiv(4)*vapor_frac+Xil(4)*(1-vapor_frac) %IM NOT SURE THIS IS ASSEMBLED CORRECTLY, BUT IM TRYING NOT TO FIXATE ON INDIVIDUAL LINES

%  BELOW THIS CONFUSES ME
%bc.dirichlet for outflux IS NOW DEFINED FOR Zi, F, Sw, and P ... (OUR PRIMARY VARIABLES)
%NOW DEFINE DIRICHLET CONDITIONS FOR INFLUX
bound_cell_in=8;
bc.influx_cells=bound_cell_in;
[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,cubic_time]=GI_flash(influxFluid,thermo,options) %DONT SEEM TO USE ANYTHING FROM THIS FLASH TEST...WHY?
for ic=1:3 %DOING WATER SEPERATELY
    bc.C_influx(ic)=influx_rate*(vapor_y(ic)*vapor_frac+liquid_x(ic)*(1-vapor_frac)); %IM NOT SURE IF THIS IS CORRECT EITHER
end
bc.water_influx=influx_rate*(liquid_x(4)*(1-vapor_frac)+vapor_y(4)*vapor_frac);

%THIS IS THE END OF SETTING UP THE CONTROLS!!!