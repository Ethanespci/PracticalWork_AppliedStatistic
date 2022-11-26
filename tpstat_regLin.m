% TP Logistic regression STAP 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
figure
load reglog_data_3.mat

plotdata(X,C)
cas = 'B'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cas 
    case 'A'
        % initialization :
        w0 = -5;
        w = randn(1,2);
        % Choice of the learning rate :
        lr = 0.01;
        % Plot the data and the line
        figure 
        subplot(2,2,1)
        plotdata(X,C,w0,w)
        title('Representation  of the initial line')

        % inference: calculate the probabilities of belonging to class 1,
        % according to the model with parameters w0 and w, for each example of X
        a = w0 + w*X; %w(1)*X(1,:) + w(2)*X(2,:) This is the same expression but less optimal
        Y = 1./(1 + exp(-a));
        %Cost function calculation
        N = length(a);
        L = - 1/N * sum(C.*log(Y + eps) + (1-C).*log(1-Y + eps)); % Add eps to avoid taking the log of 0
        % Calculation of the Cost function gradient
        dw = - 1/N * (C - Y)*X.';
        dw0 = - 1/N * sum(C - Y);
        % Update the parameter
        w = w - lr*dw;
        w0 = w0 - lr*dw0;
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning loop and monitoring: 
% We can now set up the learning loop
% The number of learning epochs is fixed by a variable. 
% The objective is to observe the evolution of certain quantities 
% during the learning process, in particular 
% the evolution of the loss function that we want to minimize. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Nepochs = 1000 %Nunber of iteration inside the loop
        Costfunc = zeros(1,Nepochs);
        Norme = zeros(1,Nepochs);
        taux_err = zeros(1,Nepochs);
        N = length(C);
        for e=1:Nepochs
            N_mplace = 0;
            a = w0 + w*X; 
            Y = 1./(1 + exp(-a)); % Calculate the probability to be in the class c=1 for each points
            L_n = - 1/N * sum(C.*log(Y + eps) + (1-C).*log(1-Y + eps)); % calculate the cost function with the current parameter of the regression
            L = L_n ;
            dw = - 1/N * (C - Y)*X.'; % calculate the gradient of the cost function
            dw0 = - 1/N * sum(C - Y); % Different calculation for w0 because the matrice X doesn't have a first line 
            w = w - lr*dw; % The parameters are changed thanks to the previous calculation of the gradient
            w0 = w0 - lr*dw0;
            Costfunc(e) = L; % This table of values keep the value of the Cost function for each iteration
            Norme(e) = sqrt(w0^2 + w*w.'); % Calculate the norm
            for k=1:N % For each dots, calculate the probability to be c=1. If it is superior to 0.5 and c=0, then the dot is on the wrong side of the delimitation line and we increase the counter
                if C(k) == 0
                    p = 1/(1 + exp(-a(k)));
                    if p > 0.5
                        N_mplace = N_mplace + 1;

                    end
                else % IT is also in the wrong position when proba<0.5 and c=1
                    p = 1/(1 + exp(-a(k)));
                    if p < 0.5
                        N_mplace = N_mplace + 1;
                    end
                end
            end
            taux_err(e) = N_mplace/N;
        end
        subplot(2,2,2)
        plotdata(X,C,w0,w) %Plot of the classification line
        title('Plot of the straight line after the regression')
        subplot(2,2,3)
        plot(Costfunc) 
        xlabel("Nbr of epoch")
        title('Evolution of the Cost function')
        
        subplot(2,2,4)
        plot(Norme)
        xlabel("Nbr of epoch")
        title('Norm of W')
% 
%         subplot(1,2,1)
%         plot(taux_err) %trace le taux d'erreur
%         xlabel("nbr d'itération")
%         title("Taux d'erreur de classification")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quadratic regression (ellipsoid shape)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'B'
        figure
        N = length(C)
        X_u = [X(1,:).^2; X(2,:).^2; X(1,:).*X(2,:); X(1,:); X(2,:)]; %Set a new vector X
        % w = randn(1,5); % Set a new vector for the parameters of the regression
        w0 = 0.5 ;
       
        subplot(2,2,1)
        plotdata(X,C,w0,w)
        title("Initial plot of the ellipse")

        Nepochs = 1000000
        lr = 0.001;
        Costfunc = zeros(1,Nepochs);
        taux_err = zeros(1,Nepochs);

        for e=1:Nepochs
            a = w0 + w*X_u ;
            Y = 1./(1 + exp(-a));
            L = - 1/N * sum(C.*log(Y + eps) + (1-C).*log(1-Y + eps));
            dw = - 1/N * (C - Y)*X_u.';
            dw0 = -1/N * sum(C - Y);
            w = w - lr*dw;
            w0 = w0 -lr*dw0;
            Costfunc(e) = L;
           
        end
        

        subplot(2,2,2)
        plotdata(X,C,w0,w)
        title("plot of the ellipse after regression")

        subplot(2,2,[3,4])
        plot(Costfunc)
        xlabel("Number of epoch")
        title('Evolution of the cost function')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for graphical representation of the data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotw(bias, w)
%%%% TO DRAW the line in case of a
%linear regression
%fplot(@(x) -(bias + w(1)*x )/w(2), [0 20]) 
%%%%% TO DRAW the ellipse in case of the quadratic regression
ezplot( sprintf('%f*x^2 + %f*y^2 + %f*x*y + %f*x + %f*y + %f = 0' ,w,bias),[-100,100,-100,100])  

end

function plotdata(MX, MY,b,w)
    neg = MY==0;
    pos = MY==1;
    plot(MX(1,neg), MX(2,neg), 'r.'); hold on;  
    plot(MX(1,pos), MX(2,pos), 'g+');
    if nargin == 4
        plotw(b,w)
    end
    xlim([0 20])
    ylim([0 20])
    hold off
end
