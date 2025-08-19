function model = model2adjacency(model)

I = model.I;

% init adjacency matrix
nnodes = I.nstates + I.npar;
adjmat = zeros(nnodes);
speciesmask = [ones(1, I.nstates), zeros(1, I.npar)];
substratemask = zeros(nnodes);

% create sym ode
X_sym = transpose(sym(I.nmstate));
par_sym = transpose(sym(I.nmpar));
dX_sym = expand(model.ode(0, X_sym, par_sym, model));

% parse sym ode
for i = 1:I.nstates
if ~isequal(dX_sym(i), sym(0))
    % read current reaction
    currstr = string(dX_sym(i));

    % identify reaction components
    min = strfind(currstr,'-');
    plu = strfind(currstr,'+');
    separators = sort([min, plu]);
    components = strings(0);
    if ~isempty(separators)
        components(end+1) = extractBetween(currstr,1,separators(1)-2);
        for sepi = 1:(length(separators) - 1)
            components(end+1) = extractBetween(currstr,separators(sepi),separators(sepi+1)-2);
        end
        components(end+1) = extractBetween(currstr,separators(end),strlength(currstr));
    else
        components = currstr;
    end

    % loop through components
    for compi = 1:length(components)
        component = components(compi);

        % identify direction of reaction
        if extract(component,1) == "+"
            isproduct = true;
            component = extractAfter(component, 2);
        elseif extract(component,1) == "-"
            isproduct = false;
            component = extractAfter(component, 2);
        else
            isproduct = true;
        end

        % identify separating operators
        mul = strfind(component,'*');
        div = strfind(component,'/');
        pow = strfind(component,'^');
        sepops = sort([mul, div, pow]);
        subcomponents = strings(0);

        % identify corresponding subcomponents
        if ~isempty(sepops)
            subcomponents(end+1) = extractBetween(component,1,sepops(1)-1);
            for subsepi = 1:(length(sepops) - 1)
                subcomponents(end+1) = extractBetween(component,sepops(subsepi)+1,sepops(subsepi+1)-1);
            end
            subcomponents(end+1) = extractBetween(component,sepops(end)+1,strlength(component));
        else
            subcomponents = component;
        end

        % label subcomponents
        species = strings(0);
        params = strings(0);
        for subcomponent = subcomponents
            if any(strcmp(I.nmstate,subcomponent))
                species(end+1) = subcomponent;
            elseif any(strcmp(I.nmpar,subcomponent))
                params(end+1) = subcomponent;
            % else
                % probably a power
                % error(['ODE component neither species nor parameter: ', char(subcomponent)])
            end
        end
        if length(params) > 1
            error("More than one parameter in reaction not supported.")
        end

        % update adjacency matrix and substrate mapping
        % contributions of other species
        if ~isempty(params) && isproduct
            for spec = species
                % update adjmat
                adjmat(I.(spec), I.nstates + I.(params)) = 1;

                % check if substrate to be able to remove arrowheads
                for substratespec = species
                    % read potential substrate reaction
                    substratestr = string(dX_sym(I.(substratespec)));
                    issubstrate = false;

                    % parse substratestr
                    % identify reaction components
                    substratemin = strfind(substratestr,'-');
                    substrateplu = strfind(substratestr,'+');
                    substrateseparators = sort([substratemin, substrateplu]);
                    substrcomponents = strings(0);
                    if ~isempty(substrateseparators)
                        substrcomponents(end+1) = extractBetween(substratestr,1,substrateseparators(1)-2);
                        for substratesepi = 1:(length(substrateseparators) - 1)
                            substrcomponents(end+1) = extractBetween(substratestr,substrateseparators(substratesepi),substrateseparators(substratesepi+1)-2);
                        end
                        substrcomponents(end+1) = extractBetween(substratestr,substrateseparators(end),strlength(substratestr));
                    else
                        substrcomponents = substratestr;
                    end
                    % loop through components
                    for substrcompi = 1:length(substrcomponents)
                        substrcomponent = substrcomponents(substrcompi);

                        % remove - and + from start
                        if extract(substrcomponent,1) == "+"
                            substrcomponent = extractAfter(substrcomponent, 2);
                        elseif extract(substrcomponent,1) == "-"
                            substrcomponent = extractAfter(substrcomponent, 2);
                        end
                
                        % identify separating operators
                        substrmul = strfind(substrcomponent,'*');
                        substrdiv = strfind(substrcomponent,'/');
                        substrpow = strfind(substrcomponent,'^');
                        substrsepops = sort([substrmul, substrdiv, substrpow]);
                        substrsubcomponents = strings(0);
                
                        % identify corresponding subcomponents
                        if ~isempty(substrsepops)
                            substrsubcomponents(end+1) = extractBetween(substrcomponent,1,substrsepops(1)-1);
                            for substrsubsepi = 1:(length(substrsepops) - 1)
                                substrsubcomponents(end+1) = extractBetween(substrcomponent,substrsepops(substrsubsepi)+1,substrsepops(substrsubsepi+1)-1);
                            end
                            substrsubcomponents(end+1) = extractBetween(substrcomponent,substrsepops(end)+1,strlength(substrcomponent));
                        else
                            substrsubcomponents = substrcomponent;
                        end
                
                        % label subcomponents
                        substrspecies = strings(0);
                        substrparams = strings(0);
                        for substrsubcomponent = substrsubcomponents
                            if any(strcmp(I.nmstate,substrsubcomponent))
                                substrspecies(end+1) = substrsubcomponent;
                            elseif any(strcmp(I.nmpar,substrsubcomponent))
                                substrparams(end+1) = substrsubcomponent;
                            % else
                                % probably a power
                                % error(['ODE component neither species nor parameter: ', char(subcomponent)])
                            end
                        end

                        % loop to substratespec components to find matches
                        if all(ismember(species, substrspecies)) && all(ismember(params, substrparams))
                            issubstrate = true;
                        end
                    end

                    % if substrate, update substratemask
                    if issubstrate
                        substratemask(I.(substratespec), I.nstates + I.(params)) = 1;
                    end
                end
            end
        end
        % species-less reactions
        if ~isempty(params)
        if isproduct
            adjmat(I.nstates + I.(params), I.(I.nmstate{i})) = 1;
        else
            adjmat(I.(I.nmstate{i}), I.nstates + I.(params)) = 1;
        end
        end
    end

end
end

model.adjmat = adjmat;
model.speciesmask = speciesmask;
model.substratemask = substratemask;

% save('adjmat.mat', 'adjmat', 'mask', 'substratemask')

end