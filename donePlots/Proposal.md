The Mystic Lake model solves an equation called `flux` using the ODE15s solver. This function is taken as an argument by the ODE solver, which feeds this function the concentration matrix (flattened into a concentration vector). The function computes the rates of reactions, diffusions, and precipitations and outputs the fluxes for each metabolite at each depth.

```matlab 
% twiddle because this is time-independent

function [conc_fluxes] = flux(~, concs_vector)

	% These are the matrices containing the concentrations and fluxes of each species
    % There are 10 species (columns) and 17 layers (1 meter depth sections)
	concs = reshape(concs_vector, [n_x, n_species]);
    conc_fluxes = zeros(n_x, n_species);

    %... some code is omitted ...

    % loop from surface to bottom
    for x = 1: n_x
        %% ... skip down to the relevant part ....

        % diffusion is based on a D_plus and a D_minus coefficients for all species
        % The following lines state that: ...
        if x > 1
            conc_fluxes(x, :) = conc_fluxes(x, :) + D_plus .* concs(x - 1, :) - D_minus .* concs(x, :);
            % this line adds what diffuses down from above and subtracts what diffuses up from this layer
        end

        if x < n_x
            conc_fluxes(x, :) = conc_fluxes(x, :) - D_plus .* concs(x, :) + D_minus .* concs(x + 1, :);
           	% this line subtracts what diffuses down from this layer, and addes what diffuses up from below
        end
    end
    % .... more code is omitted ...
end

```



