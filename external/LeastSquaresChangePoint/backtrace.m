function answer = backtrace(sq_resid, P, changepoints, N)
%************************************************************************/
%* Change Point - A program for detecting regime shifts via             */
%*                  least squares linear regression                     */
%*                                                                      */
%* Please acknowledge the program authors on any publication of         */
%* scientific results based in part on use of the program and           */
%* cite the following article in which the program was described.       */
%*                                                                      */
%* Ruggieri, E., T. Herbert, K. T. Lawrence, and C. E. Lawrence (2009), */
%* Change point method for detecting regime shifts in paleoclimatic time*/
%* series: Application to d18O time series of the Plio-Pleistocene,     */
%* Paleoceanography, 24, PA1204, doi:10.1029/2007PA001568.              */
%* Program Author: Eric Ruggieri                                        */
%*                                                                      */
%* Copyright (C) 2008   Brown University                                */
%* Brown University                                                     */
%* Providence, RI 02912                                                 */
%* Email:  Eric_Ruggieri@brown.edu                                      */
%*                                                                      */
%* This file is part of Change Point.                                   */
%*                                                                      */
%* Change Point is free software: you can redistribute it and/or modify */
%* it under the terms of the GNU General Public License as published by */
%* the Free Software Foundation, either version 3 of the License, or    */
%* (at your option) any later version.                                  */
%*                                                                      */
%* Change Point is distributed in the hope that it will be useful,      */
%* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
%* GNU General Public License for more details.                         */
%*                                                                      */
%* You should have received a copy of the GNU General Public License    */
%* along with Change Point.  If not, see <http://www.gnu.org/licenses/> */
%* or write to the Free Software Foundation, Inc., 51 Franklin Street   */
%* Fifth Floor, Boston, MA  02110-1301, USA.                            */
%************************************************************************/

% This function will back trace through the table that was created in the
% partition function and find min{f_k(v) + sq_resid(v+1,j)}
% P = Partition Function
% Indicies must match the forward step which ran from 1->v and v+1->j

last_pos=N;
answer = zeros(1,N); 
if (changepoints ~=0)                    % A homogeneous segment
    for i=changepoints:-1:2
        resid = P(i-1,1:last_pos-1)' + sq_resid(2:last_pos,last_pos);
        resid = resid - P(i,last_pos);   %The squared error you are trying to find/match
        [min_resid index] = min(resid);  %Find the position of the change point, minimum = 0 (matching values)
        answer(index+1:last_pos) = i;    %Fill in our 'prediction'
        last_pos = index;                % 1->v and v+1->j
    end
    %When i=1, the final change point
    resid= sq_resid(1,1:last_pos-1) + sq_resid(2:last_pos,last_pos)';
    [min_resid index] = min(resid);
    answer(1:index) = 0;
    answer(index+1:last_pos) =1;
end