% TP STAP 2021

clear all       
clc            
format compact  
cas = 'A';      
ALPHA = 0.05;   % Risk of error of the first kind
chemin = '';    % path for the data

switch cas
    case 'A'
        display('Height distribution (Galton)');
        noms = {'pères', 'mères', 'couples'};
        nom_fichier = sprintf('%sdataA_Tpstat.txt', chemin);
        x = load(nom_fichier);
        x(:,3) = (x(:,1) + 1.08*x(:,2))/2;
        clear figure
        for i = 1:3
            subplot(3,3,i)
            fprintf('\n*** %s ***\n', noms{i});
            tailles = x(:,i);
            % plot the histogramm
            %hist(tailles)
            [eff,xb] = hist(tailles);
            bar(xb,eff)
            hold on
            % Tracer l'histogramme cumulé, etc.
            subplot(3,3,i+3)
            effcumul = cumsum(eff);
            bar(xb,effcumul)
            title(noms{i});
            xlabel('taille en mètre');
            ylabel("nombre d'individus");
            
            variance = var(tailles); %Non biases estimator, equivalent to var(x,0)
            esperance = mean(tailles);
            ecart_type = std(tailles);
            coeffasym = skewness(tailles); %en anglais le coefficient d'asymétrie c'est skewness
            coeffaplat = kurtosis(tailles); %en anglais c'est kurtosis
            fprintf('\n espérance = %f \n variance = %f \n écart-type = %f \n coeff asymétrique = %f \n coeff aplatissement = %f \n ', esperance, variance, ecart_type, coeffasym, coeffaplat)
            
            [h, p,STAT] = chi2gof(tailles)
            subplot(3,3,i+6)
            qqplot(tailles)
        
        end
        correlation = corr(x(:,1),x(:,2))
        
    otherwise
        display('PCRq sur des données de trisomie 21');
        nom_fichier = sprintf('%sdataB.txt', chemin);
        x = load(nom_fichier);
        %the first 20 rows are TS 21+ (i.e. those who have a problem) and then the TS 21-. 
        % Average of triplicates (i.e. column averages)
        HPRTmoy = mean(x(:,1:3),2);
        TRIOBPmoy = mean(x(:,4:6),2);
        MYLIPmoy  = mean(x(:,7:9),2);
        alpha=0.05;
        % Normalise to HPRT (soustraction des ct)
        MYLIPnorm = MYLIPmoy-HPRTmoy;
        TRIOBPnorm = TRIOBPmoy-HPRTmoy;
        data = [TRIOBPnorm MYLIPnorm];
        
        % Bagplot
        groupe = [zeros(10,1); ones(10,1)];
        subplot(1,2,1)
        boxplot(TRIOBPnorm, groupe, 'labels', {'TS21+', 'TS21-'})
        title('TRIOB')
        
        subplot(1,2,2)
        boxplot(MYLIPnorm, groupe)
        
        % Analyse différentielle (comparaison des TS21+ et TS21-)
        noms = {'TRIOBP', 'MYLIP'};
        for i=1:2
            fprintf('\n*** Gène %s ***\n', noms{i});
            gene = data(:,i);
            geneplus = gene(1:10); %cardiac
            genemoins = gene(11:20); %non cardiac
            
            %normality test for small samples (n<30)
            %swtest
            
            [Hp, pValuep, Wp] = swtest(geneplus, alpha);
            [Hm, pValuem, Wm] = swtest(genemoins, alpha);
       
            subplot(2,2,i)
            qqplot(geneplus)
            subplot(2,2,i+2)
            qqplot(genemoins)
            
            %Fisher test
            if (Hp == 0 && Hm == 0)
                [hv,pv,icv,statsv] = vartest2(geneplus,genemoins)
                %Student test
                disp('** Test de Student **')
                if hv == 0 
                    [ht,pt,ict,statst] = ttest2(geneplus,genemoins)
                end
            else % Wilcoxon test
                [P,H,STATS] = ranksum(geneplus,genemoins)
                nboot = 10000
                mp = bootstrp(nboot,@mean,geneplus); 
                %bootstrap it makes nboot random drawing of our basic values 
                % and it makes the average of them and do it nboot times so 
                % we obtain a populatin of nboot =1000
                mm = bootstrp(nboot,@mean,genemoins);
                d = mp - mm;
                hist(d)
            end
            
            
        end
end
