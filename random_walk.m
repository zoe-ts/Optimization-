function [x_star,f_star] = random_walk(f,n,x,x_k)
%An implementation of Random Walk Method with Direction Exploitation for optimization 
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
e=0.05;				%minimum allowable step length 
max_iter=100;    	% maximum number of iterations
rmin=-1; 	%low r bound
rmax=1;		%upper r bound
R=2;
i=1;		%iteration number
optl=1; 	%initial step length

%Function evaluation
f_new=subs(f,x,x_k);		
f1=f_new;

%Random walk algorithm
while (optl>e) &&(i<=max_iter) 
	while R>1
		r=rmin+rand(1,n)*(rmax-rmin);	
		R=0;
		for j=1:n
			R= r(j)^2+R;
		end
	end
	R=R^(-1/2);
	u=R*r';		%unit vector
	syms l;
	g= subs(f,x, (x_k+u*l));
	newg=diff(g,l)==0;
	optl=solve(newg,l);			%optimal step length
	x_k1=x_k+(optl*u);			%new vector
    %function evaluation
	f_new=subs(f,x,x_k1);		
	if f_new<f1
		x_k=x_k1;
		f1=f_new;
		i=1;
	elseif f_new>=f1
		if i<=max_iter
			i=i+1;
		elseif i>max_iter
			l=l/2;
		end
	end 
end

x_star=x_k1;
f_star=f1;
end			






