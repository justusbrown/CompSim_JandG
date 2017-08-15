function totalFluid=addBWfluid(system, mole_fraction, m_i);
components=system.components;
Temp=system.Temp;
Cells=system.Cells;
nCell=system.nCell;

totalFluid=struct();

for i=1:nCell
    totalFluid(i).bip = zeroBIP(components);
    totalFluid(i).components = components;
    totalFluid(i).mole_fraction = mole_fraction;
    totalFluid(i).pressure=Cells.pressure(i);
    totalFluid(i).Temp=Temp;
    totalFluid(i).m_i=m_i;
end
totalFluid=totalFluid';
end

    