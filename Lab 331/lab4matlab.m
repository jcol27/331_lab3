% Plot the explained variance ratio
clear all;

% Get input
A = csvread("Eigenvalues.csv");
sum = 0;
for i = 1:size(A,1)
    sum = sum + A(i,2);
    cum(i) = sum;
end

a = size(A,1);
for i=1:size(A,1)
    x(i) = i;
    y(i) = cum(i)/sum;

end

figure(1);
scatter(x,y);
title("Plot of the Explained Variance Ratio");
xlabel("Eigenvalue index");
ylabel("Explained Variance Ratio");

% Plot the spectra on the most important principle components
C = csvread("Sorted Eigenvectors.csv");
Spectra = csvread("DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv");
Spectra = Spectra(2:57,:);
Spectra = Spectra.';
x1 = C(1,:)*Spectra(:,:);
y1 = C(2,:)*Spectra(:,:);


figure(2);
scatter(x1,y1);
title("Projection on the Two Most Important Principle Components");
xlabel("Largest");
ylabel("Second Largest");

% Plot 1500 graphs



for i=1:10
    for j=i:10
        for k=1:56
            x(i,j,k) = C(i,:)*Spectra(:,k);
            y(i,j,k) = C(j,:)*Spectra(:,k);
        end
        for k=30:56
            x(i,j,k) = C(i,:)*Spectra(:,k);
            y(i,j,k) = C(j,:)*Spectra(:,k);         
        end
    end
end


figure(3);
scatter(x(1,2,1:29), y(1,2,1:29),'r')   
hold on;
scatter(x(1,2,30:56), y(1,2,30:56),'b')
title('1st and 2nd');

figure(4);
scatter(x(2,3,1:29), y(2,3,1:29),'r')   
hold on;
scatter(x(2,3,30:56), y(2,3,30:56),'b')
title('2nd and 3rd');

figure(5);
scatter(x(3,4,1:29), y(3,4,1:29),'r')   
hold on;
scatter(x(3,4,30:56), y(3,4,30:56),'b')
title('3rd and 4th');

figure(6);
scatter(x(1,3,1:29), y(1,3,1:29),'r')   
hold on;
scatter(x(1,3,30:56), y(1,3,30:56),'b')
title('1st and 3rd');

figure(7);
scatter(x(1,4,1:29), y(1,4,1:29),'r')   
hold on;
scatter(x(1,4,30:56), y(1,4,30:56),'b')
title('1st and 4th');

figure(8);
scatter(x(2,4,1:29), y(2,4,1:29),'r')   
hold on;
scatter(x(2,4,30:56), y(2,4,30:56),'b')
title('2nd and 4th');

% Reconstruction
 Spectra = csvread("DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv");
 E = C(1:6,:);
 F = E.';
 TransformationMatrix = F*E;
 G = Spectra(10,:);
 H = Spectra(40,:);
 G = G.';
 H = H.';
 I = TransformationMatrix*G;
 J = TransformationMatrix*H;
 
 figure(9);
 plot(Spectra(1,:), G);
 hold on;
 plot(Spectra(1,:),I);
 hold off;
 title("Arabica");
 
 figure(10);
 plot(Spectra(1,:), H);
 hold on;
 plot(Spectra(1,:),J);
 hold off;
 title("Robusta");
 
 
 






