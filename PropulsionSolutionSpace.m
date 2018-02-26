alpha = linspace(0,5,6);
PRc = linspace(10,30,6);
Alt = 35000;
M0 = 1.6;
AB = 1;
S_plot = zeros(length(alpha),length(PRc));
F_plot = zeros(length(alpha),length(PRc));

for i=1:length(alpha)
    for j = 1:length(PRc)
        Sfun = @(PRf) OptiS(alpha(i),PRf,PRc(j),Alt,M0,AB);
        %[PRfopt,S] = fminbnd(Sfun,3,4,optimset('TolFun',.1,'MaxIter',10,'display','none'));
        try
        [PRfopt,S] = GoldenSection(1,5,Sfun)
        catch continue
        end
        if isreal(S)==0
            S = 0;
        end
        S_plot(i,j) = S;
        try
        F = OptiF_mdot(alpha(i),PRfopt,PRc(j),Alt,M0,AB);
        catch continue
        end
        if isreal(F)==0
            F = 0;
        end
        F_plot(i,j) = F;
    end
end
        
surf(PRc,alpha,S_plot)
xlabel('PRc')
ylabel('alpha')
zlabel('TSFC')

figure
surf(PRc,alpha,F_plot)
xlabel('PRc')
ylabel('alpha')
zlabel('F/mdot')