function estrela = StarMatrix
estrela = zeros(20);

for i = 1:20
   estrela(i,i) = 3.;
   estrela(i,20-i+1) = 1.;
end

estrela(10,1:9) = 1.;
estrela(1:9,10) = 1.;
estrela(11,12:20) = 1.;
estrela(12:20,11) = 1.;
estrela(10,10) = 11.;
estrela(11,11) = 11.;