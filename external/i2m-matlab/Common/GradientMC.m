function [h,J,time_diff]=GradientMC(m,C,hinit,Jinit,M,learning_rate,iter_thres,diff_thres)
%Gradient descent to get the h and J from m and C. hinit and Jinit are used
%as the seed. M is the number of steps in the Monte Carlo algorithm. learning_rate is the gradient step, and the descent stops
%after iter_thres steps, or if the difference between the h and J obtained
%by the model, and the empirical h and J provided as arguments, is below
%diff_thres. This difference is printed at each time step. 
N=size(m,1)

h=zeros(N,1);
J=zeros(N,N);

h=hinit;
J=Jinit;

t=0;

current_diff = 1000000;

while ((t<iter_thres)&&(current_diff>diff_thres))
    t=t+1;
    [m_th, C_th] = MC(h,J,N,M);
    current_diff=sqrt(mean((m_th-m).^2)+mean(mean((C_th-C).^2)))
    time_diff(t) = current_diff;
    if (current_diff>diff_thres)
        h = h - (learning_rate)*(m_th - m);
        J = J - (learning_rate)*(C_th - C);
    end
end

