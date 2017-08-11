function [meta, residuals] = BWgetResiduals(meta, eqs, system, linsolver_diverged)
   if nargin < 4, solver_diverged = false; end
%
    % Store the residuals for debugging and convergence testing.
    residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);

    if isempty(meta.res_history) 
        meta.res_history = zeros(system.options.nonlinear.maxIterations, numel(residuals));
    end

    meta.res_history(meta.iteration, :) = residuals;

    % Try a simple detection of oscillations, and relax the next iteration if
    % oscillations were detected.
    [oscillate stagnate] = detectNewtonOscillations(meta.res_history, system.options.cellwise, meta.iteration, system.options.nonlinear.relaxRelTol);
    if ~ linsolver_diverged,
        if oscillate
            meta.relax = max(meta.relax - system.options.nonlinear.relaxInc, system.options.nonlinear.relaxMax);
            dispif(mrstVerbose, ...
                  ['Oscillating behavior detected: Relaxation set ', ...
                   'to %.1g\n'], meta.relax);
        elseif stagnate && 0
            meta.relax = max(meta.relax + system.options.nonlinear.relaxInc, system.options.nonlinear.relaxMax);
            dispif(mrstVerbose, ...
                  ['Stagnating behavior detected: Relaxation set ', ...
                   'to %.1g\n'], meta.relax);
        end
    end
    meta.oscillate = oscillate;
    meta.stagnate  = stagnate;
    meta.linsolver_diverged = linsolver_diverged;
end
