
function [state0]=initState_JandG(rock,components,Temp,p, options,thermo);

numCells=rock.G.cells.num;
R = getR_JandG();
state.fluid.muL=1e-3;
state.fluid.muG=1e-5;
state.fluid.cl    = 4.4e-5/atm;  % Compressibility
state.fluid.p_ref = 1*atm;       % Reference pressure
totalFluid=zeros(numCells,1);

for i=1:numCells
    totalFluid(i)=addmixture(components, T, p);
    
    
[success_flag,stability_flag,Xiv,Xil,Zgas_vap,Zgas_liq,vapor_frac,cubic_time]=GI_flash(totalFluid(i),thermo,options);

totalFluid(i).Xig=Xiv(1:3); %4 components. units=MOLig/MOLg
totalFluid(i).Xio=Xil(1:3); %units=MOLio/MOLo
totalFluid(i).Xwv=Xiv(4); %units=MOLwv/MOLw
totalFluid(i).Xwl=Xil(4);
totalFluid(i).V=vapor_frac;
totalFluid(i).So=.25;
totalFluid(i).Sg=.30;
totalFluid(i).Sw=1.*ones(numCells,1)-totalFluid(i).So-totalFluid(i).Sg;
totalFluid(i).Zi=totalFluid(i).Xig.*totalFluid(i).V+totalFluid(i).Xio.*(1-totalFluid(i).V); %ALREADY Added to each cell
%totalFluid{i}.Zi=num2cell(totalFluid{i}.Zi,1);
totalFluid(i).Eo=totalFluid(i).pressure/(Zgas_liq*R*totalFluid(i).fluid.temperature); %ALREADY Added to each cell
totalFluid(i).Eg=totalFluid(i).pressure/(Zgas_vap*R*totalFluid(i).fluid.temperature); %ALREADY Added to each cell
totalFluid(i).F=(totalFluid(i).Eo.*totalFluid(i).So+totalFluid(i).Eg.*totalFluid(i).Sg);%ALREADY Added to each cell
totalFluid(i).Ew=55.5; 

%I MIGHT CHANGE EVERYTHING THAT IS A STATE TO TOTAL FLUID?
end


%CONSTANT VISCOSITIES & COMPRESSIBILITY (which I don't know if I need atm)

state0=state;

end


