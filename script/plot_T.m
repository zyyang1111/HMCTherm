close all 
clear
clc

for l = 0 : 5
    T = csvread(['temperature_layer' num2str(l) '.csv']);
    figure
    imagesc(T);
    title(['Temperature (Layer ' num2str(l) ')']);
    
    fprintf('Layer %d : maxT = %.2f\n', l, max(T(:)));
end