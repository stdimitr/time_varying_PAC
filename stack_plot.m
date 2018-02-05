function Y=stack_plot(X,step);
%
% this m-file seperate the curves with a step 
% so as to be ready to be plotted as a stack with the first one at the bottom
% Y=forplot(X,step);

Y=X;[N,p]=size(X);

for i=1:N;
  Y(i,:)= X(i,:)+ step*(i-1);
end
