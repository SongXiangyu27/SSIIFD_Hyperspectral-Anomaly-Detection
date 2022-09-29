function welcome
%WELCOME Displays DR Toolbox version information
%
%   welcome
%
% Displays DR Toolbox version information.

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7b.
% The toolbox can be obtained from http://ticc.uvt.nl/~lvdrmaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Tilburg University, 2008

    global DR_WELCOME;
    if isempty(DR_WELCOME)
        disp(' ');
        disp('   Welcome to the Matlab Toolbox for Dimensionality Reduction, version 0.7b (16-November-2008).');
        disp('      You are free to modify or redistribute this code (for non-commercial purposes), as long as a reference');
        disp('      to the original author (Laurens van der Maaten, TiCC, Tilburg University) is retained.');
        disp('      For more information, please visit http://ticc.uvt.nl/~lvdrmaaten/dr');
        disp(' ');
        DR_WELCOME = 1;
    end
