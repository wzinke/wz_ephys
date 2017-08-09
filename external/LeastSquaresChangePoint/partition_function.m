function P =partition_function(kmax,N, sq_resid)
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

% This function completes the forward recursion to calculate the partition function
% Denoted here: P(k, y_[1:j] ); in manuscript: f_k( y_[1:j] )
% The first loop is the initialization 

P=zeros(kmax,N) + Inf;

% The First Change Point %
for j=1:N               
    minimum = Inf;
    for v=1:j-1         %location of a single changepoint
        if (sq_resid(1,v) +sq_resid(v+1,j) <minimum)   
                minimum = sq_resid(1,v) +sq_resid(v+1,j);
        end
    end
    P(1,j) = minimum;
end

% Now loop on successively more changepoints %
for kk=2:kmax;
    for j=kk+1:N;  
        minimum = Inf;
        for v=kk:j-1;    %location of final changepoint
            if(P(kk-1,v) +sq_resid(v+1,j) <minimum)
                minimum = P(kk-1,v) +sq_resid(v+1,j);
            end
        end                     %end of v loop
        P(kk, j) = minimum;
    end                         %end of j loop 
end                             %end of kk loop
