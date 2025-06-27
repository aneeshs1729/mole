% SPDX-License-Identifier: GPL-3.0-only
% 
% Copyright 2008-2024 San Diego State University Research Foundation (SDSURF).
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% LICENSE file or on the web GNU General Public License 
% <https://www.gnu.org/licenses/> for more details.
%
% ------------------------------------------------------------------------
    % Returns a (m+2) x (m+1) one-dimensional mimetic divergence operator

    % Input validation
    function D = div(k, m, dx)
% DIV  Mimetic 1-D divergence operator using Björck-Pereyra Vandermonde solver
%    D = DIV(k,m,dx) returns the (m+2)×(m+1) sparse divergence matrix of
%    order-k accuracy on a uniform grid of spacing dx.
%
%    Compared with the original version, this code
%    ──────────────────────────────────────────────
%      • avoids building full Vandermonde matrices;
%      • relies on the O(k²) Björck-Pereyra algorithm (primal form);
%      • is numerically stabler for large k.

    %% sanity checks ------------------------------------------------------
    assert(k >= 2 && mod(k,2)==0, 'k must be an even integer ≥ 2');
    assert(m >= 2*k+1, 'Need at least 2k+1 cells for order-k scheme');

    nRows = m + 2;         % = #faces
    nCols = m + 1;         % = #cells
    D     = sparse(nRows, nCols);

    %% central rows (interior stencil) ------------------------------------
    % stencil abscissae relative to the cell face
    s = (1/2 - k/2) + (0:k-1);             % length-k vector
    d=zeros([1,k]);
    d(2)=1;
    coeffs = vandSolveBP(s, d);      % 1st-derivative moments
     j = 1;
    for i = k/2+1 : nRows - k/2
        D(i, j:j+k-1) = coeffs;
        j = j + 1;
    end
    %% boundary closures --------------------------------------------------
    p = k/2 - 1;      % #extra rows to build at each end
    q = k + 1;        % stencil width for boundaries
    d=[d 0];
    for i = 1:p 
            s = (1/2 - i) + (0:q-1);       % shifted stencil
            c = vandSolveBP(s, d);   % derivative coefficients
                  % left block
            D(1+i, 1:q) = c;
            c=-c(end:-1:1);
            D(nRows-i, nCols-q+1:nCols)=c;
    end
    D = D / dx;   % scale by grid spacing
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = vandSolveBP(alpha, b)
% vandSolveBP  Björck-Pereyra O(n²) solver for V(alpha)*a = b.
%
%   alpha : vector of abscissae (length n)
%   b     : right-hand side  (length n)
%   a     : solution         (length n)
%
% The implementation follows the two-sweep triangular scheme published in
% Björck & Pereyra (1970).

    n = numel(alpha) - 1;          % index runs 0..n
    a = b(:).';                    % work in row-vector form

    %% forward sweep  (Step (i) – triangular L factor)
    for k = 1:n
        for j = n+1:-1:k+1
            j;
            a(j) = a(j) - alpha(k)*a(j-1);
        end
    end

    %% backward sweep (Step (ii) – triangular U factor)
    for k = n:-1:1
        for j = k+1:n+1
            a(j) = a(j) / (alpha(j) - alpha(j-k));
        end
        for j = k:n
            a(j) = a(j) - a(j+1);
        end
    end
end