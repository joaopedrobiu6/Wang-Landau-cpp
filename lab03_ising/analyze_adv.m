% script to analyze MC data in matlab
% Martin Lind�n 20030923
% Petter Johansson 2015
graphics_toolkit("gnuplot")

% Set input data files to read from
datafile={'data4','data6','data8','data10','data12','data16','data20'};
mark={'.b-','*r-','^k-','dg-','.m-','*b-','^r-','.k:','.g:','*m-'};

% The critical temperature
Tc = 2/log(1+sqrt(2));

% Set the critical exponents here
nu = 1;
beta = 1/8;

% Plot data from input files into a single figure
figure(1); clf; hold on;
for k = 1:length(datafile);
   % Load data
   try
       d = load(datafile{k});
   catch
       printf("Could not read file '%s'.\n", datafile{k});
       continue;
   end

   % Some parameters
   T = d(:, 1);
   m = d(:, 2);
   mabs = d(:, 3);
   g = d(:, 7);
   L = d(:, 8);

   xi = L.^(beta/nu).*mabs; % Magnetization scaled by L and beta
   eta = L.^(1/nu).*(T-Tc); % Temperature scaled by L and nu

   plot(T, m, mark{k});
end

d = load('exact_solution.data');
k = k+1;
T_exact = d(:, 1);
e_exact = d(:, 2);
c_exact = d(:, 3);
mabs_exact = d(:, 4);
% Uncomment to add something from the exact solution
%plot(T_exact, mabs_exact, mark{k});

xlabel('T');
ylabel('m'); % change
title('Ising model, short run');

% Uncomment to set axis limits manually
%axis([0 4 0 1]);

% Uncomment to print figure
%print('fig.eps','-deps');

% Uncomment to save figure data as a text file
%xy=[x y];
%save('fig.plot','-ascii','xy');

