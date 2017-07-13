components=addComponents({'CH4','C4H10','H2O'});
Nc=single(numel(components));;
Np=single(3);

rho_c=ones(Nc, Np); %mass density comp i phase j
CH4=1;C4H10=2;H2O=3;
E_c=ones(Nc, Np); %Molar density comp i phase j
E=ones(Np,1); %Molar Denisty phase j
oil=1;gas=2;water=3;
x=ones(Nc,Np);%Molar fraction comp i phase j
mass=ones(Nc,Np)
pv=1; %Temporary placeholder for pore volume

for j=1:Np
    for i=1:Nc
        mass=components(i).MW*x(i,j);
        rho_c(i,j)=mass/pv;
    end
end



for j=1:Np
    for i=1:Nc
        E_c(i,j)=rho_c(i,j)/components(i).MW; %molar density_comp related to MW and mass density_comp
        E(j)=sum(E_c(i,j)); %j=o,g i=1...Nc %%Molar Density of phase
        x(i,j)=E_c(i,j)/E(j); %mole fraction of component i in phase j
    end
end
%{
m=ones(i,j);
m_setup1=ones(i,j);
m_setup2=ones(j);

for j=1:Np
for i=1:Nc
m_setup1(i,j)=components(i).MW*x(i,j);
end
end

for j=1:Np
    for i=1:Nc
        m_setup2(j) += components(i).MW*x(i,oil) %WHY IS THIS GIVING TROUBLE?
    end
end

 for j=1:Np
     for i=1:Nc
         m=m_setup1(i,j)/m_setup2(j)
     end
 end
 %}

pot=ones(Np,1);
P=ones(Np,1);
rho_p=ones(Np,1);
%SETUP potential---------------------------------------------
g=9.8; %temporary placeholder
h=1; %temporary placeholder
for j=1:Np
    pot(j)=P(j)-rho_p(j)*g*h;
end
%total mass variable F-------------------------------------
S=ones(Np,1);
F=E(oil)*S(oil)+E(gas)*S(gas);

%mass fractions O&G---------------------- note L+V=1
L=E(oil)*S(oil)/F 
V=E(gas)*S(gas)/F

%total mol fraction----------- note sum(zi)=1 and xioEoSo+xigEgSg=Fzi
for i=1:Nc
z(i)=L*x(i,oil)+(1-L)*x(i,gas);
end

%Transmitibilities--------------------------------
k=1; %temporary placeholder
kr=ones(Np,1);
mu=ones(Np,1);

for j=1:Np
T(j)=E(j)*kr(j)/mu(j)*k;
end

for j=1:Np
    for i=1:Nc
        T_c(i,j)=x(i,j)*E(j)*kr(j)/mu(j)*k;
    end
end%USE MRST

%cons. mass components--------------------------------------------------------
div=@(x) x ;%placeholder
q=ones(Np,1);
phi=0.3 %placeholder
dt=1%placeholder

for i=1:Nc
massCons_c(i)=phi*F*z(i)/dt-div(T_c(i,oil)*div(pot(oil))+T_c(i,gas)*div(pot(gas)))-x(i,oil)*q(oil)-x(i,gas)*q(gas);
end

%Cpns mass TOTAL---------------------
massCons_=phi*F/dt-div(T(oil)*div(pot(oil))+T(gas)*div(pot(gas)))-q(oil)-q(gas)


%cons mass water---------------------
massCons_T=phi*E(water)*S(water)/dt-div(T(water)*div(pot(water)))-q(water)


%Saturation State Equations---------------------------------
satEq=F*(L/E(oil)+(1-L)/E(gas))+S(water)-1

%%THIS DIFF SYSTEM CONSISTS OF 2nC+2 EQ'S FOR THE SAME AMOUNT OF PRIMARY
%%UNKNOWNS XIO,l,f,sW,P, I=1...nC-1








