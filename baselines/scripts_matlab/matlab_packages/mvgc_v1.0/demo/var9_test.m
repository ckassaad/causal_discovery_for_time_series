%% var9_test
%
% Create VAR coefficients for 9-node test network
%
% <matlab:open('var9_test.m') code>
%
%% Syntax
%
%     A = var9_test
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _output_
%
%     A          VAR coefficients matrix for 9-node test network
%
%% Description
%
% Returns coefficients |A| for a VAR(3) based on a 9-node test network with causal
% architecture:
%
% <<var9test.png>>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function A = var9_test

n = 9;
p = 3;

A = zeros(n,n,p);

A(:,:,1) = [
    0.0114         0    0.4088    0.4236         0         0         0         0         0;
   -0.3394   -0.0192         0    0.2799         0         0         0    0.3085         0;
         0         0    0.0194    0.3437         0         0         0         0         0;
         0    0.4302         0    0.0188         0         0         0         0         0;
         0         0    0.2851         0    0.0027         0    0.3160         0         0;
         0         0         0         0    0.3039   -0.0020         0         0         0;
         0         0         0         0         0    0.3030   -0.0186         0         0;
         0         0         0         0         0         0         0   -0.0084    0.3477;
         0         0         0         0         0         0    0.3037    0.2810   -0.0208
];

A(:,:,2) = [
    0.0148         0    0.2590    0.1965         0         0         0         0         0;
    0.2232   -0.0070         0    0.1779         0         0         0   -0.2383         0;
         0         0   -0.0058    0.2008         0         0         0         0         0;
         0    0.2103         0    0.0012         0         0         0         0         0;
         0         0    0.1597         0    0.0065         0    0.1989         0         0;
         0         0         0         0    0.2062    0.0177         0         0         0;
         0         0         0         0         0    0.1895   -0.0008         0         0;
         0         0         0         0         0         0         0   -0.0032    0.1381;
         0         0         0         0         0         0    0.1947   -0.1718    0.0068
];

A(:,:,3) = [
    0.0076         0    0.1434    0.1787         0         0         0         0         0;
    0.1082   -0.0065         0    0.1351         0         0         0    0.2021         0;
         0         0    0.0050    0.1826         0         0         0         0         0;
         0    0.2090         0   -0.0037         0         0         0         0         0;
         0         0    0.1943         0    0.0097         0   -0.0698         0         0;
         0         0         0         0    0.1347   -0.0081         0         0         0;
         0         0         0         0         0    0.2270   -0.0004         0         0;
         0         0         0         0         0         0         0   -0.0014    0.1397;
         0         0         0         0         0         0    0.2020    0.1417   -0.0022
];
