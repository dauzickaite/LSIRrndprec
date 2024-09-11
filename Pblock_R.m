function  Px = Pblock_R(x,convprc,R,m,u,upr)
%%%%
% Preconditioner P product with x, where
% P*x = |  I            0        | * |x1| ,
%       |  0          R^(-1)     |   |x2|  
% 
% I is m x m identity and R is n x n upper triangular.
% Computations are performed in precision upr and rounded to lower precision
% u if convprc is set to true.
%%%%

switch upr 
    case 'single'
        x = single(x);
        R = single(R);
    case 'double'
        x = double(x);
        R = double(R);
    case 'quad'
        x = mp(x,34);
        R = mp(R,34);
end


x2 =  R\x(m+1:end,:);

if convprc
    switch u 
        case 'single'
            Px = single([x(1:m,:);x2]);
        case 'double'
            Px = double([x(1:m,:);x2]);
    end
else
    Px = [x(1:m,:);x2];
end
