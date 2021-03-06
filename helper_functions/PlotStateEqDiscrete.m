function PlotStateEqDiscrete(InVar, OutVar, Cain_vect, Fc_vect, Az, Bz, Cz, Dz, titleText, plotStyle)
    global Ca T CAin Fc colour_vect;
    colour_vect = [0.905, 0.094, 0.698;
                0.278, 0.819, 0.137;
                0.905, 0.772, 0.094;
                0.109, 0.152, 0.949;
                0.258, 0.960, 0.847; ];
            
    vect = [Cain_vect;Fc_vect];

    for iter = 1:1:length(vect)

%         if (strcmp(InVar, 'CAin') == 1)
%             CAin = vect(iter);
%         elseif (strcmp(InVar, 'Fc') == 1)
%             Fc = vect(iter);
%         end

        if (strcmp(OutVar, 'Ca') == 1)
            y_plot = 1;
        elseif (strcmp(OutVar, 'T') == 1)
            y_plot = 2;
        end

%         [y, t] = rk4(dCaIn, dTIn, Ca, T, step);
        n = 1e3;
        x = [CAin;Fc].*ones(2, n);
        y = zeros(2, n);
        for k = 2:1:n
            x(:, k) = Az*x(:, k-1)+Bz*[Cain_vect(iter);Fc_vect(iter)];
            y(:, k) = Cz*x(:, k) + Dz*[Cain_vect(iter);Fc_vect(iter)];
        end
            
        plot((1:n), y(y_plot, :), plotStyle, 'Color', colour_vect(iter, :), 'LineWidth', 1.5);
        title(titleText)
        xlabel('t [min]')

        if (strcmp(OutVar, 'Ca') == 1)
            ylabel('Ca [kmol/m^3]');
        elseif (strcmp(OutVar, 'T') == 1)
            ylabel('T [K]');
        end

        hold on
    end
end

