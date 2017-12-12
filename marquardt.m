function [x_star,f_star] = marquardt(f,n,x,x_k)
%An implementation of Marquardt's Method for optimization 
%
%Inputs: 
%   f: objective function
%   n: number of dimensions
%   x: symbolic vector of n size
%   x_k: initial point
%Outputs:
%   x_star: optimal point
%   f_star: optimal value of objective function
%
%Author: Zoe Ts. 2017 


% Initialization
e=10^-2;  		% tolerance
k=1;			% iteration
a=10^4;  		% tolerance
c1= 1/4;
c2= 2;
I= eye(n);		% identity matrix 

%Gradient of the objective function
Grad_f = sym(zeros(n,1));
for i = 1:n
   Grad_f(i) = gradient(f, x(i,1)); 			
end

%Hessian matrix
H_f = sym(zeros(n, n));
for i = 1:n
   for j = 1:n
       H_f(i, j) =  diff(Grad_f(i), x(j));		
   end
end

Grad_fk=subs(Grad_f,x,x_k);
df_k=norm(Grad_fk,Inf);
tol=abs(df_k);

% Marquardt's algorithm
while (tol>=e)
    Grad_fk=subs(Grad_f,x,x_k);
	dH=H_f+a.*I;						% inverse of the hessian matrix
    s=-dH^(-1)*Grad_fk;					%search direction
    f_k=subs(f,x,x_k);
	Grad_fk=subs(Grad_f,x,x_k);
	tol=abs(Grad_fk);
	x_k1=x_k+s;							%new estimation
	f_k1=subs(f,x,x_k1);
	while f_k1>=f_k
		a=c2*a;
		x_k1=x_k+s;						%new estimation if f_k1>=f_k
		f_k1=subs(f,x,x_k1);
	end
	if f_k1<f_k 
		a=c1*a;
		k=k+1;							%new iteration if f_k1<f_k
	end
	x_k=x_k1;
end

x_star=x_k;
f_star=subs(f,x,x_k);
end



