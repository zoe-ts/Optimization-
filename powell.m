function [x_star,f_star] = powell(f,n,x,x_k)
%An implementation of Powell's Method for optimization 
%
%Inputs: 
%   f: objective function
%   n: number of dimensions
%   x: symbolic vector of n size
%   x_k: initial point
%Outputs:flet
%   x_star: optimal point
%   f_star: optimal value of objective function
%
%Author: Zoe Ts. 2017 

% Initialization
e=0.01;					%probe length
i=1;					%iteration
f_new=1;
fplus=0;
fminus=0;

% Initial guess
xvalue= zeros(n,n+1);	%stores values of x_k	
s=[0;1];				%initial search direction

jj=1;
%Powell's algorithm
while (fplus<f_new) || (fminus<f_new)

% Cycle1: Univariate algorithm
while i<=3
	if i<n-1
		j=i+1;
		s=zeros(n,1);
		s(j,1)=1;					
	elseif i>=n
		s=zeros(n,1);
		s(jj,1)=1;
		jj=jj+1;
		if jj>n
		jj=1;
		end
	end
	f_new=subs(f,x,x_k);				%calculating f
	x_kplus=x_k+e*s;
	fplus=subs(f,x,x_kplus);			%calculating f+	
	x_kminus=x_k-e*s;
	fminus=subs(f,x,x_kminus);			%calculating f-	
	%correct direction s and step 
	syms l ;
	if fplus<f_new
		g= subs(f, x, (x_k+s*l));
		newg=diff(g,l)==0;
		optl=solve(newg,l);		%optimal step length
		x_k1=x_k+optl*s;
		xvalue(1:n,i)=x_k1;
		x_k=x_k1;
	elseif fminus<f_new
		g= subs(f, x, (x_k-s*l));
		newg=diff(g,l)==0;
		optl=solve(newg,l);		%optimal step length
		x_k1=x_k-optl*s;
		xvalue(1:n,i)=x_k1;
			x_k=x_k1;
	end
	i=i+1;
end

% Cycle2: Pattern search		
	s_p=x_k-xvalue(1:n,1);              %pattern direction
	f_new=subs(f,x,x_k);				%calculating new f
	x_kplus=x_k+e*s_p;
	fplus=subs(f,x,x_kplus);			%calculating f+	
	x_kminus=x_k-e*s_p;
	fminus=subs(f,x,x_kminus);			%calculating f-	
	%correct pattern direction s and step 
	syms l;
	if fplus<f_new
		g= subs(f, x, (x_k+s_p*l));
		newg=diff(g,l)==0;
		optl=solve(newg,l);		%optimal step length
		x_k1=x_k+optl*s_p;
		f_new=subs(f,x,x_k1);	%calculating new f
	else
		g= subs(f, x, (x_k-s_p*l));
		newg=diff(g,l)==0;
		optl=solve(newg,l);		%optimal step length
		x_k1=x_k-optl*s_p;
		f_new=subs(f,x,x_k1);   %calculating new f
	end

end %end Powell's algorithm

x_star=x_k1;
f_star=subs(f,x,x_star);
end