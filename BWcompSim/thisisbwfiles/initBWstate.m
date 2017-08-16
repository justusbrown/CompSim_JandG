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
muL=1e-3;
muG=1e-5;
cl    = system.cl;  % Compressibility
p_ref = system.p_ref;       % Reference pressure

%THE USER DEFINES THE INITIAL ZI, NOT MI? JB 8/15
for i=1:nCell
    totalFluid(i).Zi=totalFluid(i).mole_fraction;
end

Xig=[]; %4 componentss. units=MOLig/MOLg
Xio=[]; %units=MOLio/MOLo
V=[];
So=[];
Sg=[];
Sw=[];
Eo=[]; 
Eg=[]; 
Ew=[]; 
rhoL=[];
rhoG=[];
m_w=[];
m_i=[];


for i=1:nCell    
    %{
    CAN WE DO THIS INSTEAD? THIS MAY BE RETARDED BUT JUST A THOUGHT
    m_i=[m_i; totalFluid(i).m_i]; %IS THIS USER DEFINED IN ADD BW FLUID? DO
    WE NEED BOTH THIS AND ZI?
    state.m_i=m_i;
    state.m_i=num2cell(state.m_i,1);
    totalFluid(i).Zi=totalFluid(i).m_i./sum(totalFluid(i).m_i);
    %}
    
    [success_flag,stability_flag,Xiv,Xil, vapor_frac, Zgas_vap,Zgas_liq,cubic_time]=GI_flash(totalFluid(i), thermo,options);

    m_i=[m_i; totalFluid(i).m_i]; %IS THIS USER DEFINED IN ADD BW FLUID? DO WE NEED BOTH THIS AND ZI (CAN WE DO THE ABOVE INSTEAD?)
    state.m_i=m_i;
    state.m_i=num2cell(state.m_i,1);

    
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
    
    totalFluid(i).Eo=totalFluid(i).pressure/(Zgas_liq*R*totalFluid(i).Temp); %ALREADY Added to each cell
    Eo=[Eo;totalFluid(i).Eo];
    state.Eo=Eo;

    totalFluid(i).Eg=totalFluid(i).pressure/(Zgas_vap*R*totalFluid(i).Temp); %ALREADY Added to each cell
    Eg=[Eg;totalFluid(i).Eg];
    state.Eg=Eg;
    
    totalFluid(i).Ew=55.5/system.litre; 
    Ew=[Ew;totalFluid(i).Ew];
    state.Ew=Ew;

    totalFluid(i).m_w=totalFluid(i).Ew*(1 - sum(totalFluid(i).m_i)*(1-totalFluid(i).V)/totalFluid(i).Eo - sum(totalFluid(i).m_i)*totalFluid(i).V/totalFluid(i).Eg);
    m_w=[m_w; totalFluid(i).m_w];
    state.m_w=m_w;
    
    totalFluid(i).rhoL=totalFluid(i).Eo*sum(MW'.*totalFluid(i).Xio);
    rhoL=[rhoL; totalFluid(i).rhoL];
    state.rhoL=rhoL;

    totalFluid(i).rhoG=totalFluid(i).Eg*sum(MW'.*totalFluid(i).Xig);
    rhoG=[rhoG; totalFluid(i).rhoG];
    state.rhoG=rhoG;
    
end

state.totalFluid=totalFluid;

for i=1:nCell
   state.p(i)=totalFluid(i).pressure;
end
state.p=state.p';

state0=state;

end






    
    