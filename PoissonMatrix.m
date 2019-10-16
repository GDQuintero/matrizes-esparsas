function poisson = PoissonMatrix(n,grid)

poisson = 5*eye(n);

for i = 1:n-1
   poisson(i+1,i) = -1;
    poisson(i,i+1) = -1;
end

for i = 1:n-grid
   poisson(grid+i,i) = -1; 
   poisson(i,grid+i) = -1;
end



