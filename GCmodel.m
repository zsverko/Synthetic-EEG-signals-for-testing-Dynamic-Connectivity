function [GC, A1, A2, A12, e1, e2, e12] = GCmodel(data, order)
%GCmodel performs the computation of Granger Caucality for two row vectors provided in the data field.
% (c) Peter Rogelj: peter.rogelj@upr.si
% The function estimates Granger causality for predicting the first row of
% data using additional second row of data.
% Inputs:
%   data - matrix with two rows each presenting one variabel. Prediction is made for the first one
%   order - the number of history time points used for prediction in VAR models
% Outputs:
%   GC - the resulting Granger Causality
%   A  - n parameters for the univariate VAR model for the first raw of the data (n=order)
%   AB - 2*n parameters for the bivariate VAR model for predicting the first row using both rows of data
%   e1 - estimated univariate prediction error
%   e2 - estimated bivariate prediction error

X = zeros(2*order, size(data,2)-order); % contains shifted data for uni and bivariate VAR
for n=1:order  
    X(n,:) = data(1,1+order-n:end-n); %univariate +
    X(order+n,:) = data(2,1+order-n:end-n); %bivariate +
end
x1=data(1,1+order:end); 
x2=data(2,1+order:end); 
x12=data(1:2,1+order:end); 

% bivariate VAR model
%A12 = x12/X;
A12 = x12*pinv(X);
e12 = x12 - A12*X;
% univariate VAR 1
%A1 = x1/X(1:n,:); 
A1 = x1*pinv(X(1:n,:)); 
e1 = x1 - A1*X(1:n,:);
% univariate VAR 2
%A2 = x2/X(order+(1:n),:); 
A2 = x2*pinv(X(order+(1:n),:)); 
e2 = x2 - A2*X(order+(1:n),:);
% Granger causality
GC=log( [var(e1) var(e2)]./[ var(e12(1,:)) var(e12(2,:)) ] );

%disp(A);
%disp(AB);

end

