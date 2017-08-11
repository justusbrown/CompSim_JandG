%C is    C = C1_interior + C1_exterior - (C2_interior + C2_exterior);
%gr 07/20

%grad_JandG function

function dval = BWgrad(val, bc_val, nf, C, ...
                     is_dirichlet_faces1, dirichlet_faces1, ...
                     is_dirichlet_faces2, dirichlet_faces2)

   signed_bc_val = sparse(nf, 1);
   signed_bc_val(dirichlet_faces1) = - bc_val(is_dirichlet_faces1);
   signed_bc_val(dirichlet_faces2) = + bc_val(is_dirichlet_faces2);
   dval = C*val + signed_bc_val;

end
