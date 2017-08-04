 
function [ops]=setupBWsystem(rock, bc);

   cf = rock.G.cells.faces(:,1);
   nf = rock.G.faces.num;
   nc = rock.G.cells.num;
   
   %%
   % Compute the half, and then the full, transmissibilities.
   T=rock.T;
   T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
   
   %%
   %Setup the discrete divergence operator.
   [C1,C2,C_div,div]=BWdivOp(rock.G, cf, nf, nc);
   
   %%
   %Setup the discrete gradient operator.
   %We compute the differences of cell values
   % across each face. It is a linear mapping from cells' to faces' values.
   N = double(rock.G.faces.neighbors);
   index = 1:nf;
   interior = prod(N, 2)~=0;

   C1_interior = sparse(index(interior), N(interior, 1), ones(nnz(interior), 1), nf, nc);
   C2_interior = sparse(index(interior), N(interior, 2), ones(nnz(interior), 1), nf, nc);

   %%
   % Compute the boundary contribution to the gradient operator. They corresponds to the
   % external faces where Dirichlet conditions are given. We are careful to use the
   % correct signs.
   is_dirichlet_faces1 = N(bc.dirichlet.faces, 1) ~= 0;
   is_dirichlet_faces2 = ~is_dirichlet_faces1;

   dirichlet_faces1 = bc.dirichlet.faces(is_dirichlet_faces1);
   dirichlet_faces2 = bc.dirichlet.faces(is_dirichlet_faces2);

   C1_exterior = sparse(index(dirichlet_faces1), ...
                        N(dirichlet_faces1, 1), ...
                        ones(numel(dirichlet_faces1), 1), nf, nc);
   C2_exterior = sparse(index(dirichlet_faces2), ...
                        N(dirichlet_faces2, 2), ...
                        ones(numel(dirichlet_faces2), 1), nf, nc);
                    
   %%
   % The gradient operator is the sum of internal and boundary contributions.
   C = C1_interior + C1_exterior - (C2_interior + C2_exterior);

   pressure_bc = sparse(nf, 1);
   pressure_bc(dirichlet_faces1) = - bc.dirichlet.pressure(is_dirichlet_faces1);
   pressure_bc(dirichlet_faces2) = + bc.dirichlet.pressure(is_dirichlet_faces2);

   p_grad = @(p)(C*p + pressure_bc);  

   grad = @(val, bc_val)(grad_JandG(val, bc_val, nf, C, ...
                                is_dirichlet_faces1, dirichlet_faces1, ... 
                                is_dirichlet_faces2, dirichlet_faces2));

   %%
   % Set up the gravity term.
   z = rock.G.cells.centroids(:, 3);
   fz = rock.G.faces.centroids(:, 3);
   dz = grad(z, fz);
   
   %% 
   % Define function to compute concentrations on faces using cells' values and upwind
   % directions. Here, the _nature_ of the boundary faces (Dirichlet or note) are not
   % allowed to change, although the _values_ on these faces may change.
   faceConcentrations = @(flag, conc_c, bc_conc) ...
       BWfaceConcentrations(flag, conc_c, bc_conc, N, interior, dirichlet_faces2, ...
                          dirichlet_faces1,  bc, nf, nc);
%%
%ops will be the operator struct containing everything from setupBWsystem.
ops.div=div; ops.p_grad=p_grad; ops.grad=grad; ops.dz=dz;
ops.faceConcentrations=faceConcentrations;
ops.N=N; ops.G=G;
end