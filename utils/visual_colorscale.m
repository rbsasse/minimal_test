function c = visual_colorscale(n, varargin)
%VISUAL_COLORSCALE Gradient color scale equispaced in hsv
%
% Use as follows to generate colours:
% colorscale(n)
% colorscale(n, 'hue', [min max])
% colorscale(n, 'saturation', saturation)
% colorscale(n, 'value', value)
%
% PARAMETERS:
% - n
% - hue in [0 1]x[0 1]
% - range (default [0.1 0.9])
% - saturation [0 1] (default 0.5)
% - value in [0 1] (default 0.8)
%
% RETURNS: nx3 rgb matrix
%
% DEPENDS: none
%
% FAMILY: user_level, utility, graphics
%
% EXAMPLES:
% n = 10;
% cols = colorscale(n, 'hue', [0.1 0.8], 'saturation' , 1, 'value', 0.5);
%
%for aa = 1:10;
%     plot(1:10, (1:10) + aa, 'Color', cols(aa,:), 'Linewidth',2)
%     hold on
%end;
%
% % plot a matrix
% v = transpose(1:10);
% set(gca, 'ColorOrder', colorscale(5));
% set(gca,'NextPlot','replacechildren')
% plot(v, [v, v+1, v+2, v+ 3, v+4, v+5]) ;
% % set as default for all plots
% set(groot,'defaultAxesColorOrder', visual_colorscale(10))
%

% convoluted logic to parse optional arguments
p = inputParser;
p.addRequired('n', @isnumeric);
p.addOptional('hue', [0.1 0.9], @(x) length(x) == 2 & min(x) >=0 & max(x) <= 1);
p.addOptional('saturation', 0.5, @(x) length(x) == 1);
p.addOptional('value', 0.8, @(x) length(x) == 1);

% parse what was given
p.parse(n, varargin{:});

c = hsv2rgb([transpose(linspace(p.Results.hue(1), p.Results.hue(2), p.Results.n)), ...
    repmat(p.Results.saturation, p.Results.n, 1), repmat(p.Results.value, n,1) ]);

end
