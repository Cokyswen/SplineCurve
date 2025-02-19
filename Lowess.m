## Copyright (C) 2025
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} Lowess (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <cxw@cxw-pc>
## Created: 2025-02-19
function y_smooth = Lowess(data, alpha, order)
    n = length(data);
    y_smooth = zeros(n, 1);
    for i = 1:n
        % 计算每个点的权重
        distances = abs(data(:, 1) - data(i, 1));
        max_distance = max(distances);
        weights = (1 - (distances / (alpha * max_distance)).^3).^3;
        weights(distances > alpha * max_distance) = 0;

        % 加权最小二乘拟合
        X = [ones(n, 1), data(:, 1) - data(i, 1)];
        if order > 1
            for p = 2:order
                X = [X, (data(:, 1) - data(i, 1)).^p];
            end
        end
        W = diag(weights);
        beta = (X' * W * X) \ (X' * W * data(:, 2));
        y_smooth(i) = beta(1);  % 拟合值
    end
endfunction
