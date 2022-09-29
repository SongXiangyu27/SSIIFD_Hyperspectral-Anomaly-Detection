function reconX = run_data_through_autoenc(network, X)
%RUN_DATA_THROUGH_AUTOENC Summary of this function goes here
%   Detailed explanation goes here

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7b.
% The toolbox can be obtained from http://ticc.uvt.nl/~lvdrmaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Tilburg University, 2008

    % Initialize some variables
    n = size(X, 1);
    no_layers = length(network);
    middle_layer = ceil(no_layers / 2);

    % Run data through autoencoder
    activations = [X ones(n, 1)];
    for i=1:no_layers
        if i ~= middle_layer
            activations = [1 ./ (1 + exp(-(activations * [network{i}.W; network{i}.bias_upW]))) ones(n, 1)];
        else
            activations = [activations * [network{i}.W; network{i}.bias_upW] ones(n, 1)];
        end
    end
    reconX = activations(:,1:end-1);
