   %% 
   %%IMPORTANT NOTE: gr 07/20
   %I THINK WE MIGHT NEED TWO OF THESE BECAUSE WE WILL NEED TWO DIFFERENT
   %FLUXES, ONE WITH FZi and one WITH JUST F.
   %ACTUALLY, Im not sure we will after all
   %
   % Define function to compute concentrations on faces using cells' values and upwind
   % directions. Here, the _nature_ of the boundary faces (Dirichlet or note) are not
   % allowed to change, although the _values_ on these faces may change.
   %
   function conc_f = faceConcentrations_JandG(flag, conc_c, bc_conc, N, interior, ...
                                                 dirichlet_faces2, dirichlet_faces1, ...
                                                 bc, nf, nc)
   index        = (1:nf)';
   upCell       = N(:, 2);
   upCell(flag) = N(flag, 1);
   
      % On the interior cell we use upwind
   Mint = sparse(index(interior), upCell(interior), 1, nf, nc);

   logical_dirichlet_faces1 = zeros(nf, 1);
   logical_dirichlet_faces1(dirichlet_faces1) = 1;
   logical_dirichlet_faces1 = logical(logical_dirichlet_faces1);
   logical_dirichlet_faces2 = zeros(nf, 1);
   logical_dirichlet_faces2(dirichlet_faces2) = 1;
   logical_dirichlet_faces2 = logical(logical_dirichlet_faces2);

   external_faces1 = N(:,2)==0;
   external_faces2 = N(:,1)==0;


   % On the exterior faces where no Dirichlet conditions are given we take the value given in
   % the interior cell.
   Mext1 = sparse(index(external_faces1 & ~logical_dirichlet_faces1), ...
                  N(external_faces1 & ~ logical_dirichlet_faces1, 1), 1, nf, nc);
   Mext2 = sparse(index(external_faces2 & ~logical_dirichlet_faces2), ...
                  N(external_faces2 & ~ logical_dirichlet_faces2, 2), 1, nf, nc);

   % On the Dirichlet boundary cells we use upwind, taking the values from boundary
   % conditions when needed We assume flag is logical

   assert(islogical(flag), 'Upstream indices must be given as logical');
   Mdir1 = sparse(index(flag & logical_dirichlet_faces1), ...
                  N((flag & logical_dirichlet_faces1), 1), 1, nf, nc);
   Mdir2 = sparse(index(~flag & logical_dirichlet_faces2), ...
                  N((~flag & logical_dirichlet_faces2), 2), 1, nf, nc);

   M = Mint + Mext1 + Mext2 + Mdir1 + Mdir2;

   % Values of saturation from Dirichlet boundary conditions.
   dconc_all = sparse(nf, 1);
   dconc_all(bc.dirichlet.faces) = bc_conc; 
   dconc = sparse(nf, 1);
   dconc(flag & logical_dirichlet_faces2)  = dconc_all( flag & logical_dirichlet_faces2);
   dconc(~flag & logical_dirichlet_faces1) = dconc_all(~flag & logical_dirichlet_faces1);
   
   conc_f = M*conc_c + dconc;
   
end
%BEGINNING OF NONLINEAR SOLVER PARAMETERS, DIDN'T CHANGE, DOESN'T SEEM TO
%NEED TO BE


   
  