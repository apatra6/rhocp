clc; clear;
%% MATLAB script to plot IPF

CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.9 2.9 2.9], 'mineral', 'Ferrite','color', 'red'),...
  crystalSymmetry('m-3m', [3.6 3.6 3.6], 'mineral', 'Austenite', 'color', 'green')};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% original input microstructure
fname = 'example_2_out.ang';

ebsd = EBSD.load(fname,CS,'interface','ang','convertEuler2SpatialReferenceFrame','setting 1');
[grains,ebsd.grainId] = calcGrains(ebsd);
plot(grains); % phase map
% exportgraphics(gcf,'phase_map_example_2.png','Resolution',300);

ipfKey = ipfColorKey(grains('Ferrite'));
ipfKey.inversePoleFigureDirection = zvector;
color = ipfKey.orientation2color(grains('Ferrite').meanOrientation);
figure(2)
plot(grains('Ferrite'),color);
exportgraphics(gcf,'ipf_map_example_2.png','Resolution',300);
figure(3)
plot(ipfHSVKey(ebsd('Ferrite').CS))
exportgraphics(gcf,'ipf_map_example_2_key.png','Resolution',300);
