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
## @deftypefn {} {@var{retval} =} CatmullRomSpline (@var{alpha}, @var{numSamplesPerSpline}, @var{points})
##
## @seealso{}
## @end deftypefn

## Author:  <cxw@cxw-pc>
## Created: 2025-02-18

function C = CatmullRomSpline(alpha, samplesPerNum, points)
    pointCount = size(points, 1);
    pointDimensions = size(points, 2);

    if pointCount < 4
        error('There must be at least 4 control points!');
    end

    if pointCount > 4
        splinesNum = pointCount-3;

        C = zeros(samplesPerNum*splinesNum, pointDimensions);
        for i=1:splinesNum
            spline = CatmullRomSpline(alpha, samplesPerNum, points(i:i+3, :));

            stopIndex = i * samplesPerNum;
            startIndex = stopIndex - samplesPerNum + 1;

            C(startIndex:stopIndex, :) = spline;
        end
    end

    if pointCount == 4
        P0 = points(1, :);
        P1 = points(2, :);
        P2 = points(3, :);
        P3 = points(4, :);

        t0 = 0;
        t1 = getT(t0, alpha, P0, P1);
        t2 = getT(t1, alpha, P1, P2);
        t3 = getT(t2, alpha, P2, P3);

        t = linspace(t1, t2, samplesPerNum)';

        A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1;
        A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2;
        A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3;

        B1 = (t2-t)/(t2-t0).*A1 + (t-t0)/(t2-t0).*A2;
        B2 = (t3-t)/(t3-t1).*A2 + (t-t1)/(t3-t1).*A3;

        C  = (t2-t)/(t2-t1).*B1 + (t-t1)/(t2-t1).*B2;
    end
endfunction
