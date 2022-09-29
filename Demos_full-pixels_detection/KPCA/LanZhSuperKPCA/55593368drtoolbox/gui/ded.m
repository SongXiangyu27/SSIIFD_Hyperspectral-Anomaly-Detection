function ded(cc,hs,he)
% cc - edit number
% hs - handle to plot
% he - handle to edit

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7b.
% The toolbox can be obtained from http://ticc.uvt.nl/~lvdrmaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Tilburg University, 2008

data=get(hs,'UserData');
ld=length(data(1,:));

nd=str2num(get(he,'string'));
if nd>ld
    nd=ld;
end
if nd<1
    nd=1;
end
nd=round(nd);

set(he,'string',num2str(nd));
drawnow;

switch cc
    case 1
        set(hs,'Xdata',data(:,nd));
    case 2
        set(hs,'Ydata',data(:,nd));
    case 3
        set(hs,'Zdata',data(:,nd));
end

drawnow;