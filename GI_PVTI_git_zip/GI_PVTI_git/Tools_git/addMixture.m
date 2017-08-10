function mixture =  addMixture(components, T_K, p_Pa)
% This function creates a mixture structure

mixture = struct();
mixture.bip = zeroBIP(components);
mixture.components = components;
n = length(components);
mixture.Zi = 1./(n*ones(1,n));
mixture.pressure = p_Pa; %[Pa]
mixture.Temp = T_K; % [K]