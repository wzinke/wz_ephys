function changepoint(handles)
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

%
% Outline of the code is as follows:
% 1) Load the data from the GUI and then the input file
% 2) Remove long-term trends if desired
% 3) Build the Model
% 4) Compute Squared Error for each possible substring: sq_resid(i,j)
% 5) Build the Partition Function: P(k,j) [f_k(j) in manuscript]
% 6) Decide upon the number of change points and then Backtrace to find optimal solution
% 7) Display the results
%
% Description of Parameters, input, and output can be found in ReadMe file
%

warning off all;
%*****************
%1) Load the Parameters from the GUI 
filename=(get(handles.filename_input, 'String'));
num_chgpts = str2num(get(handles.num_chgpts_input,'String'));
kmax = str2num(get(handles.kmax_input,'String'));
dist = str2num(get(handles.dist_input,'String'));

arg=get(handles.detrend_checkbox,'Value');
if (arg)
    detrend='d';
else detrend ='';
end
arg=get(handles.skip_checkbox,'Value');
if (arg)
    skip='s';
else skip ='n';
end

J=get(handles.weight_buttons,'SelectedObject');
switch get(J, 'Tag')   % Get Tag (or String) of selected object
    case 'argU'
      %execute this code when unweighted radiobutton is selected
      weight='u';
    case 'argW'
      %execute this code when weighted radiobutton is selected
      weight='w';
    otherwise
       % Code for when there is no match.
end

J=get(handles.regressor_buttons,'SelectedObject');
switch get(J, 'Tag')   % Get Tag (or String) of selected object
    case 'file_radiobutton'
      %execute this code when file radiobutton is selected
     create='f';
     omegas=[];
    case 'create_radiobutton'
      %execute this code when create radiobutton is selected
      create='c';
      omegas=[str2num(get(handles.w1, 'String')) str2num(get(handles.w2, 'String')) ...
          str2num(get(handles.w3, 'String')) str2num(get(handles.w4, 'String')) str2num(get(handles.w5, 'String'))];
      omegas=sort(omegas); 
      pause(0.01);
    case 'poly_radiobutton'
        create='p';
        
    otherwise
       % Code for when there is no match.
end

rem=filename;
while true
    [fname rem] = strtok(rem,'\');
    if isempty(rem)
        break;
    end
end
fname = strtok(fname,'.');
if (skip=='s')
    fout=[fname num2str(num_chgpts) 'chgpts_' skip '.txt'];
else 
    fout =[fname num2str(num_chgpts) 'chgpts_' detrend weight create '.txt'];
end

out=fopen(fout, 'w');
fprintf(out, ['Results of the Least Squares Change Point Algorithm on ' fname '.txt using ' num2str(num_chgpts) ' change points\n']);
fprintf(out, '\nCaveats: ');
fprintf(out, '\n1) In the least squares setting, there is no rigorous way to determine the optimal number of change points');
fprintf(out, '\n2) There is also no way to characterize the uncertainty surrounding the number of change points');
fprintf(out, '\n3) The optimal solution shown is the best in terms of squared residual error.');
fprintf(out, '\n4) There is no way to characterize the uncertainty surrounding the optimal solution.  ')
fprintf(out, '\n\tIt is possible that the true timing of the change may not be the exact time found, but relatively close');
fprintf(out, ['\n5) Keep in mind that this solution is one of N^' num2str(num_chgpts) ' possible solutions']);
fprintf(out, '\n6) Given the large number of possible solutions, it is very likely that no one solution is significantly better than any other solution');
fprintf(out, '\n7) To obtain your model, multiply the coefficients by the regressors in the indicated intervals and then add');

if (skip == 'n')

    fprintf(out, ['\n\nMaximum number of change points considered: ' num2str(kmax)]);
    fprintf(out, ['\nMinimum distance between change points: ' num2str(dist)]);
    fprintf(out, ['\nNumber of change points in output: ' num2str(num_chgpts)]);
%**************
% 1) The d18O data files used were in the form: column 1 = time, 
% column 2 = proxy value, column 3 through M = regressors

vars=dlmread(filename, '', 1,0);   %Labels in first row
[N, M]=size(vars);  %the number of data points and regressors 

if(~any(vars(:,M)))
    M=M-1;
end
%This takes care of extra column of zeros at end

t=vars(:,1);
y=vars(:,2);

%**********
% 2) Detrend? - An exponential detrend was used on the d18O data sets

if (detrend=='d')
% A Linear Model with non-polynomial terms.  This method used a 'design
% matrix' and was found on: mathworks.com
data=vars(:,2);
X = [ones(size(t))  exp(-t/1000)  t.*exp(-t/1000)];
a=X\data;
y = vars(:,2) - (a(1) + a(2)*exp(-t/1000) +a(3)*t.*exp(-t/1000));
fprintf(out, '\nTrend removed from data via exponential function');
clear data X;
set(handles.status, 'String', 'Data Detrended... Running');
pause(0.01);
end


% ************
% 3) Build the Model: The original model used three sinusoids of different 
% frequencies, plus a constat. Each column is one predictor, each row is 
% the predictor set at a given time point

if ( create =='c')
    fprintf(out, '\nSinusoids were created at Periodicities of ')
    X=[];
    for i=1:length(omegas)
        X=[X cos(2*pi/omegas(i)*t(:)) sin(2*pi/omegas(i)*t(:))];
        fprintf(out, [num2str(omegas(i)) ' ']);
    end
    X=[X ones(N,1)];
    fprintf(out, '.  \nRegression Coefficients are for cosine then sine at each periodicity shown above.');
elseif (create =='p')
    order =str2num(get(handles.poly_order, 'String'));
    fprintf(out, ['\nA polynomial model of order ' num2str(order) ' was created']);
    X=[];
    for i=order:-1:0
        X=[X t(:).^i];
    end
else
    X=vars(:,3:M);
end
%************
%4) Calculte the Squared Error

sq_resid =zeros(N,N) +Inf;

if (weight=='w')
    fprintf(out, '\nWeighted Regression was used');
    for i=1:N
        if (mod(i,100) ==0)
            set(handles.status, 'String', ['Running...' num2str(i) ' Completed']);
            pause(0.01);
        end
        for j=i+1:N
            if(t(j)-t(i)>dist)  %this accounts for non-equally spaced data points
                hat_beta=inv(X(i:j,:)'*X(i:j,:))*X(i:j,:)'*y(i:j);
                sq_resid(i,j)= (y(i:j)-X(i:j,:)*hat_beta)'*(y(i:j)-X(i:j,:)*hat_beta);
                % For weighted least sqares only 
                sq_resid(i,j) = sq_resid(i,j)/var(y(i:j));
                %
            end
        end
    end
else
    fprintf(out, '\nOrdinary (Unweighted) Linear Regression was used');
    for i=1:N
        if (mod(i,100) ==0)
            set(handles.status, 'String', ['Running...' num2str(i)]);
            pause(0.01);
        end
        for j=i+1:N
            if(t(j)-t(i)>dist)  %this accounts for non-equally spaced data points
                hat_beta=inv(X(i:j,:)'*X(i:j,:))*X(i:j,:)'*y(i:j);
                sq_resid(i,j)= (y(i:j)-X(i:j,:)*hat_beta)'*(y(i:j)-X(i:j,:)*hat_beta);
            end
        end
    end
end
set(handles.status, 'String', 'Squared Error Calculated');
pause(0.01);
%Note 1: There are likely more efficient ways to compute this, but it depends
%on the number and type of predctors. 
%Note 2: The above double loop is also the place where you can consider 
%multiple models. Create X1, X2, X3, etc. and then
%calculate sq_resid1, sq_resid2, sq_resid3, etc.  Next, find 
%min{sq_resid1, sq_resid2, sq_resid3, etc}  over all possible models. 
%Just make sure to keep track of which of your models has the minimal 
%squared error, or simply recalculate in the backtrace step

%***********
% 5) The Partition Function
P =partition_function(kmax,N, sq_resid);
set(handles.status, 'String', 'Partition Function Built');
pause(0.01);

%*****************
% 6) Decide upon the number of change points: There is no rigorous way to
% decide on the 'optimal' number of change points in the least squares
% setting.  Instead, look for the 'break in the curve' plotted below to
% decide how many you want.

else
    load changepoint_result
    fprintf(out, '\n\nA previous output of Change Point was loaded for this analysis');
    fprintf(out, ['\nMaximum number of change points considered: ' num2str(kmax)]);
    fprintf(out, ['\nMinimum distance between change points: ' num2str(dist)]);
    fprintf(out, ['\nNumber of change points in output: ' num2str(num_chgpts)]);
end %of the skip parameter

homogen = sq_resid(1,N); % no change points
sq_sum = [homogen P(:,N)'];
sq_sumA = zeros(1,kmax);  %the change in squared residual with the addition  
        %of each successive change point
for i=1:kmax
    sq_sumA(i) = sq_sum(i)-sq_sum(i+1);
end

axes(handles.axes1); bar(sq_sumA); 
title('Reduction in Squared Error with Each Additional Change Point');
ylabel('Change in Squared Residual');
xlabel('Number of Change Points');

%Remember, bar graph includes possibility of homogeneous segment

%****************
%6) Backtrace: Given a number of change points, what is the optimal
%solution?

answer = backtrace(sq_resid, P, num_chgpts, N);
%answer will be a vector indicating segments, changepoint at boundary
% ex: 000000000111111111222222222222233333333333333333444444444444.... etc

%****************
% 7) Display the results 

xx=zeros(1,N) + mean(y)-4*std(y);  %indicates change point locations in final plot
myplot=zeros(N,1);  %plot of final model
start=1;
parameters =[];
t(:) = round(t(:)); %in case time points are not integers

for i=1:N-1
    
    if (answer(i)~=answer(i+1))
        hat_beta=inv(X(start:i,:)'*X(start:i,:))*X(start:i,:)'*y(start:i);
        myplot(start:i) = X(start:i,:)*hat_beta;
        %calculate R^2 values for this interval
        resid = (y(start:i)-myplot(start:i))'*(y(start:i)-myplot(start:i));
        sq_error = (y(start:i)-mean(y(start:i)))'*(y(start:i)-mean(y(start:i)));
        
        R2 = 1-resid/sq_error;
        parameters = [parameters; t(start) t(i) hat_beta' R2];
        
        start=i+1;
        xx(i) = xx(i) +.25; 
        
    end
end

% The final sub-interval
hat_beta=inv(X(start:N,:)'*X(start:N,:))*X(start:N,:)'*y(start:N);
myplot(start:N) = X(start:N,:)*hat_beta;

%calculate R^2 values for this interval
resid = (y(start:N)-myplot(start:N))'*(y(start:N)-myplot(start:N));
sq_error = (y(start:N)-mean(y(start:N)))'*(y(start:N)-mean(y(start:N)));

R2 = 1-resid/sq_error;
parameters = [parameters; t(start) t(N) hat_beta' R2];

%calculate R^2 for the entire data set, if desired
resid = (y(1:N)-myplot(1:N))'*(y(1:N)-myplot(1:N));
sq_error = (y(1:N)-mean(y(1:N)))'*(y(1:N)-mean(y(1:N)));

total_R2 = 1-resid/sq_error;

figure(1); plot(t, y, t, myplot, t, xx)
title('your title - data, model and change point locations');
xlabel('xlabel'); ylabel('ylabel');

if(skip=='n')  %new results, so save them!
   save changepoint_result sq_resid P parameters t y X kmax dist N;
end

set(handles.status, 'String', ['Program Complete. Total R^2 is: ' num2str(total_R2)]);

fprintf(out, ['\nTotal R^2 = ' num2str(total_R2) '\n']);
    
fprintf(out, '\nInterval \t\t Regression Coefficients');
for i=1:length(hat_beta) %The number of regressors
    fprintf(out, '\t');
end
fprintf(out, 'R^2 for Interval\n');

dlmwrite(fout, parameters, 'delimiter', '\t', 'precision', 4, '-append');

fclose('all');

end
