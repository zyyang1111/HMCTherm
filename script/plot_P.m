close all 
clear
clc

for l = 0 : 5
    if l < 4
        P = csvread(['power_mem_layer' num2str(l) '.csv']);
        figure
        imagesc(P);
        title(['Mem Layer ' num2str(l)]); 
        fprintf('Mem Layer %d: totalP = %.2e\n', l, sum(P(:)));
    elseif l < 5
        P = csvread(['power_logic_layer.csv']);
        figure
        imagesc(P);
        title(['Logic Layer']); 
        fprintf('Logic Layer: totalP = %.2e\n', sum(P(:)));
    else
        P = csvread(['power_processor_layer.csv']);
        figure
        imagesc(P);
        title(['Processor Layer']); 
        fprintf('Processor Layer: totalP = %.2e\n', sum(P(:)));
    end
end