classdef Regression < handle
    % A curve fit to a set of data points, based on regression with
    % Gaussian (default) or sigmoid (optional) basis functions. 
    
    properties
        x; 
        t; 
        phi; % basis functions
        w; % weights
    end
    
    methods 
        function r = Regression(x, t, centres, width, lambda, varargin)
            % x: samples of an independent variable 
            % t: corresponding samples of a dependent variable
            % centres: a vector of function centres (should have similar
            %   range of values as x)
            % width: width of the functions (e.g. sigma of Gaussians) 
            % lambda: regularization parameter
            % useSigmoids (optional): if 1, use sigmoid instead of Gaussian
            %   basis
            
            r.x = x;
            r.t = t;
            
            if ~isempty(varargin) && varargin{1} > 0
                r.phi = makeSigmoidFunctions(centres, width);
            else 
                r.phi = makeGaussianFunctions(centres, width);
            end
            
            r.w = findOptimalWeights(x, t, r.phi, lambda);
        end
        
        function y = eval(r, x)
            % x: a new (or multiple samples) of the independent variable 
            % y: the value of the curve at x
            
            y = zeros(size(x));
            for i = 1:length(r.phi)
                y = y + r.w(i) * r.phi{i}(x);
            end
        end
    end
end

function phi = makeGaussianFunctions(centres, width)
    phi = cell(length(centres), 1);
    for i = 1:length(centres)
        phi{i} = @(x) exp(-(x-centres(i)).^2 / (2*width^2));
    end
end

function phi = makeSigmoidFunctions(centres, width)
    phi = cell(length(centres), 1);
    for i = 1:length(centres)
        phi{i} = @(x) 1./(1+exp(-2*(x-centres(i))/width));
    end    
end

function w = findOptimalWeights(x, t, phi, lambda)
    Phi = zeros(length(x), length(phi));
    for i = 1:length(x)
        for j = 1:length(phi)
            Phi(i,j) = phi{j}(x(i));
        end
    end
    w = (lambda*eye(length(phi)) + Phi' * Phi) \ (Phi' * t);
end
