% NanobubbleDigitizer data analysis package licence details:
% 
% Copyright Soumyadeep Paul (soumyadeep.paul@thml.t.u-tokyo.ac.jp)
% 
% This file is part of NanobubbleDigitizer data analysis package.
% 
% NanobubbleDigitizer is a MATLAB package designed to detect discrete 
% bubble events from transient current traces. The bubble features
% extracted from the current data using this code can be used to analyze
% nanoscopic bubble dynamics confined within solid-state nanopores. This
% package is released under the GNU GPL.
% 
% The copyright of arb is held by Soumyadeep Paul.
% 
% NanobubbleDigitizer is released under the GNU GPL.  NanobubbleDigitizer 
% is free software: you can redistribute it and/or modify it under the 
% terms of the GNU General Public License (version 3) as published 
% by the Free Software Foundation.You should have received a copy of the 
% GNU General Public Licence along with arb (see file licence.txt after 
% unpacking).  If not, see <http://www.gnu.org/licences/>.
% 
% NanobubbleDigitizer is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
% Public Licence for more details.
% 
% For full details of arb's licence see the licence file in the main 
% directory.


%transition probability matrix

%blockage_duration_array



tw = 1000*waiting_time_1Darray;
tb = 1000*blockage_duration_1Darray;
NC = dipCurrent_1Darray;
CC = riseCurrent_1Darray;


%correlation coefficient between 4 features (tw, tb, NC, CC)
A = [tw tb NC CC];
R = corrcoef(A);


figure
heatmap(R, 'FontSize', 18);
colormap default
ax = gca;
ax.XData = ["Waiting time" "Blockage duration" "Nucleation current" "Collapse current"];
ax.YData = ["Waiting time" "Blockage duration" "Nucleation current" "Collapse current"];
ax.FontSize = 18;