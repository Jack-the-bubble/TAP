function plotLegend(InVar, vect)
    colour_vect = [0.905, 0.094, 0.698;
                0.278, 0.819, 0.137;
                0.905, 0.772, 0.094;
                0.109, 0.152, 0.949;
                0.258, 0.960, 0.847; ];
    h = zeros(size(colour_vect, 1));
    h(1) = plot(NaN, NaN, '-', 'Color', colour_vect(1, :));
    h(2) = plot(NaN, NaN, '-', 'Color', colour_vect(2, :));
    h(3) = plot(NaN, NaN, '-', 'Color', colour_vect(3, :));
    h(4) = plot(NaN, NaN, '-', 'Color', colour_vect(4, :));
    h(5) = plot(NaN, NaN, '-', 'Color', colour_vect(5, :));
    equals = ' = ';
    legend(char(strcat(InVar, equals, {' '}, num2str(vect(1)))), char(strcat(InVar, equals, {' '}, num2str(vect(2)))), ...
        char(strcat(InVar, equals, {' '}, num2str(vect(3)))), char(strcat(InVar, equals, {' '}, num2str(vect(4)))), ...
        char(strcat(InVar, equals, {' '}, num2str(vect(5)))));
end
