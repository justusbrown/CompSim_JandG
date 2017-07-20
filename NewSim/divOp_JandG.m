%DIVERGENCE OPERATOR

function [C1,C2,C,div]=divOp_JandG(G, cf, nf, nc);


   % Set up the discrete divergence operator, |div|. It sums up all signed faces'
   % values in each cell.
   % 
   
   N = double(G.faces.neighbors);
   index = 1:nf';
   faces1 = N(:, 1) ~= 0;
   faces2 = N(:, 2) ~= 0;
   C1  = sparse(index(faces1), N(faces1, 1), ones(nnz(faces1),1), nf, nc);
   C2  = sparse(index(faces2), N(faces2, 2), ones(nnz(faces2),1), nf, nc);
   C = C1 - C2;
   div  = @(x)C'*x;
   
end


