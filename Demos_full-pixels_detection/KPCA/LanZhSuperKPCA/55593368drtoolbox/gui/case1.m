% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7b.
% The toolbox can be obtained from http://ticc.uvt.nl/~lvdrmaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Tilburg University, 2008

% if now parameters: parameters make invisible:   

% handles to make invisible:
hnds={'na','nat','uipanel_tp','prp','prpt','kpd','kpdt','kpR','kpRt','kgs','kgst','uipanel_k','uipanel_ft','sig','sigt','uipanel_ei','prc','prct','k','kt','mi','mit','wl','t','tt','ka'};

for hc=1:length(hnds)
    eval(['set(handles.' hnds{hc} ',''visible'',''off'')']);
end

set(handles.sm,'Enable','on');

update_kernel_uipanel;

update_type_uipanel;