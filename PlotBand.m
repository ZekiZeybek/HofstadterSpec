
function PlotBand()
for q = 2:25
    for p = 1:q
        if gcd(p,q) == 1
            PlotOne(p,q,2)
        end
    end
end

function PlotOne(p,q,lambda)

[H0] = H_0(p,q,lambda);  H0_E = sort([eig(H0)]);
[Hpi] = H_pi(p,q,lambda); Hpi_E = sort([eig(Hpi)]);

%[H0_E';Hpi_E']
%B = repmat(1/2,1,5)

line(repmat(p/q,1,q), [H0_E';Hpi_E'], "color", "b");

function [H0] = H_0(p,q,lambda)

    
if q == 1
    H0 = [2 + lambda];
elseif q == 2
    H0 = [-lambda 2; 2 lambda]; 
else
    n = 1:q;
    D = diag(lambda*cos(2*pi*n*p/q));
    H0 = D + diag(ones(1,q-1),1) + diag(ones(1,q-1),-1);
    H0(q,1) = 1; H0(1,q) = 1;
end    


function [Hpi] = H_pi(p,q,lambda)

if q == 1
    Hpi = [-lambda-2];
elseif q == 2
    Hpi = [0 0 ; 0 0];
else
    n = 1:q;
    D = diag(lambda*cos(2*pi*n*p/q - pi/q));
    Hpi = D + diag(ones(1,q-1),1) + diag(ones(1,q-1),-1);
    Hpi(q,1) = -1; Hpi(1,q) = -1;
end
    
