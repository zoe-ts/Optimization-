function [x_star,f_star] = quasi_bfgs(f,n,x,x_k)
%An implementation of Quasi Newton BFGS Method for optimization 
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
e= 0.01;  		% tolerance
tol= 10;		% initial tolerance value
i= 1;			% iteration

%Gradient of the objective function
Grad_f = sym(zeros(n,1));
for i = 1:n
   Grad_f(i) = gradient(f, x(i,1)); 			
end
%Identity matrix (initial estimate as an inverse of Hessian)
B= eye(n);		

% BFGS algorithm
while tol>e
    Grad_fk=subs(Grad_f,x,x_k);
    %search direction
    s=-B*Grad_fk;
    %optimal step length
	syms l;
	g= subs(f,x, (x_k+s*l));
	newg=diff(g,l)==0;
	optl=solve(newg,l);	
    %new estimation
    x_k1=x_k+optl*s;					

	Grad_fk=subs(Grad_f,x,x_k);
  
	Grad_fk1=subs(Grad_f,x,x_k1);
    df_k1=norm(Grad_fk1,Inf);
	tol=abs(df_k1);
    if tol<e
		x_k=x_k1;
		break;
    end
    %updated B matrix
	d=x_k1-x_k;
	g=Grad_fk1-Grad_fk; 
	B = B+(1+((g'*B*g)./(d'*g)))*((d*d')./(d'*g))-(d*g'*B)./(d'*g)- (B*g*d')./(d'*g);	
    x_k=x_k1;
    i=i+1;
end

x_star=x_k;
f_star=subs(f,x,x_k);
end