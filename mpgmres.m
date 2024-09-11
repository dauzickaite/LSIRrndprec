function [x,resnorm,relres,iter,flag] = ...
    mpgmres(A,b,x0,tol,maxit,restart,Pright,Pleft,u,uA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves Ax=b or the preconditioned system Pleft*A*Pright*y = Pleft*b, 
% where Pright*y=x via FGMRES in 4 precision:
% - uA for products with A (single/double/quad);
% - uR for products with Pright; (preset in Pright)
% - uL for products with Pleft; (preset in Pleft)
% - u for all other computations (single/double).
% Pright (Pleft) is given as a function handle such that Pright(x,true)
% returns Pright^(-1)*x computed in uR and converted to u and Pright(x,false)
% returns Pright^(-1)*x computed and stored in uR.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mp.Digits(64);


n = length(b);

switch u 
    case 'single'
        A = single(A);
        b = single(b);
        x0 = single(x0);
        
    case 'double'
        A = double(A);
        b = double(b);
        x0 = double(x0);
end

Axprod = @(x,convprc) A_matvec(x,convprc,A,u,uA);
Pl =Pleft;
Pr = Pright;

    
x0 = Pr(x0,false); % right prec
r = Pl(b,true) - Pl(Axprod(x0,false),true);

switch u 
    case 'single'
        
        r = single(r);
        beta0 = single(norm(r));

        x = single(x0);

        relres = single(zeros(maxit*restart+1,1));
        iter = single(zeros(maxit,1));   

        
    case 'double'
        r = double(r);
        beta0 = double(norm(r));

        x = double(x0);

        relres = double(zeros(maxit*restart+1,1));
        iter = double(zeros(maxit,1));                  
end

resnorm(1)=norm(r);

relres(1) = beta0/beta0;
tot_it = 0;


for out_it = 1:maxit
    flag = 1;
    
    switch u 
        case 'single'

            Z = single(zeros(n,restart));
            V = single(zeros(n,restart+1));
            H = single(zeros(restart+1,restart));

            c = single(zeros(restart,1));
            s = single(zeros(restart,1));
            
            g = single(zeros(restart,1));

        case 'double'

            Z = double(zeros(n,restart));
            V = double(zeros(n,restart+1));
            H = double(zeros(restart+1,restart));

            c = double(zeros(restart,1));
            s = double(zeros(restart,1));
            
            g = double(zeros(restart,1));

    end

    
    % make sure that r has the right precision
    r = Pl(b,true) - Pl(Axprod(x,false),true);
        
    beta = norm(r);
    V(:,1) = r/beta;
    g(1) = beta;


    for i=1:restart
        
        tot_it = tot_it +1;

        z = Pr(V(:,i),false);
        w = Pl(Axprod(z,false),true);
        
                
        switch u 
            case 'single'
                Z(:,i) = single(z);
            case 'double'
                Z(:,i) = double(z);
        end
             
        for j = 1:i
            H(j,i) = V(:,j)'*w;
            w = w - H(j,i)*V(:,j);
        end

        H(i+1,i) = norm(w);
        V(:,i+1) = w/H(i+1,i);


        % solve y = argmin|| beta*e1 - H*y || by applying Given's rotations
        % apply rotations to the new column
        for j=1:i-1
            H(j:j+1,i) = [c(j) s(j); -s(j) c(j)]*H(j:j+1,i);
        end

        [c(i),s(i)] = givens(H(i,i),H(i+1,i));
        
        H(i,i) = [c(i) s(i)]*H(i:i+1,i);
        H(i+1,i) = 0.0;
        
        g(i:i+1) = [c(i); -s(i)]*g(i);
      
        relres((out_it-1)*restart+i+1) = abs(g(i+1))/beta0; 
        %resnorm((out_it-1)*restart+i+1) = norm(b-A*x);
        
        % if tolerance is satisfied - update x_k 
        if relres((out_it-1)*restart+i+1) <= tol
            y = H(1:i,1:i)\g(1:i);
            x = x + Z(:,1:i)*y;
            %x = x + Pr(V(:,1:i)*y,true);
            flag = 0;
            iter(out_it) = i;
            break
        end 
        iter(out_it)  = i;

    end
     
    
    if  relres((out_it-1)*restart+i+1) <= tol
        if norm(b -  A*x)<= tol*(norm(b)+norm(full(A))*norm(x)) %checking real residual before restarting
            %norm(mp(b) -  mp(A)*mp(x))<= tol*(norm(mp(b))+norm(mp(full(A)))*norm(mp(x))) %checking real residual before restarting
            break
        end 
    else
        y = H(1:i,1:i)\g(1:i);
        x = x + Z(:,1:i)*y;
        %x = x + Pr(V(:,1:i)*y,true);
    end
    
end

relres = relres(1:tot_it+1);



function Ax = A_matvec(x,convprc,A,u,uA)
% (vct) <- [vct,log,mtx,str,str]
% Compute a matrix-vector product with the A in precision uA,
% and round to precision u if convprc is set to true.

switch uA 
    case 'single'
        x = single(x);
        A = single(A); 
        
    case 'double'
        x = double(x);
        A = double(A); 
        
    case 'quad'
        mp.Digits(34);
        x = mp(x);
        A = mp(A); 
end

if convprc
    switch u 
        case 'single'
            Ax = single(A*x);

        case 'double'
            Ax = double(A*x);
    end
else
    Ax = A*x;
end

function [c,s] = givens(a,b)
if b == 0
    c=1;
    s=0;
elseif abs(b)>abs(a)
    tau = a/b;
    s=1/((1+tau^2)^0.5);
    c=s*tau;
else
    tau = b/a;
    c = 1/((1+tau^2)^0.5);
    s = c*tau;
end
   






    