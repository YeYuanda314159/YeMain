function continuation_Emin(nelx, nely, volfrac, sd, bc, g, fileID, save_to_folder, same, sameEmin)
% nelx: number of elements on x axis
% nely: number of elements on y axis
% volfrac: volume fraction of material to total area
% sd1: filter size for thresholding
% sd2: filter size for solving equation
% sd3: filter size for the perimeter term
% bc:  takes values in {'cantilever_rb', 'cantilever_central', 'mbb'},
%      enforces BCs  (boundary condition)
% g: coefficient of the perimeter term
% continuation: 1 indicates to use x as the initial guess
%               0 indicates to use constant density as the initial guess
% x: initial guess, if continuaton is 0, can be anything.
% fileID: the opened file to log outputs
% If you don't want to save the results, just let fileID to be negative 
% and the results will be displayed in the command window
% sameEmin: calculate the compliance each iteration using this specified value as Emin
    Emins =  [0.5, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001, 0.0001];
    allenergies=[];
    allenergiessame = [];
    %% start continuation
    tic; %开始计时
    sds = [sd];
    if same %whether calculate compliance with respect to the same Emin
        [y, loop, c, x, energies, energiessame] = topthr(nelx, nely, volfrac, Emins(1), g, sd, bc, 0, 1, fileID, sameEmin);
        %function [y, loop, c,  x, energies, energiessame]=topthr(nelx, nely, volfrac, frac, g, sd, bc, continuation, x, fileID, sameEmin)
        allenergiessame(1:length(energiessame)) = energiessame;
    else
        [y, loop, c, x, energies, energiessame] = topthr(nelx, nely, volfrac, Emins(1), g, sd, bc, 0, 1, fileID);
    end
    try
        saveas(gcf, strcat(mydir, meshstr, string(0), '.png'));
    catch err
        disp(err.message);
        save_to_folder = false;
    end
    allenergies(1:loop) = energies;
    numofiterations = loop;
    %%
    for i = 2:length(Emins) %依次计算每个人造材料对应的优化结果
        if same
            [y, loop, c, x, energies, energiessame]=topthr(nelx, nely, volfrac, Emins(i), g, sd, bc, 1, x, fileID, sameEmin);
            allenergiessame(numofiterations + 1 : numofiterations + length(energiessame)) = energiessame;
        else
            [y, loop, c, x, energies, energiessame]=topthr(nelx, nely, volfrac, Emins(i), g, sd, bc, 1, x, fileID);
        end
        allenergies(numofiterations + 1 : numofiterations + loop) = energies; 
        %记录柔度能的变化
        numofiterations = numofiterations + loop;
        if save_to_folder
            saveas(gcf, strcat(mydir, meshstr, string(i-1), '.png'));
        end
        figure;
        plot(numofiterations + 1 : numofiterations + loop, energies, 'o-'); axis tight;
    end  
    toc;    %结束计时
    figure; imshow(1-x); caxis([0,1]); axis equal; axis off; set(gca,'Units','normalized','Position',[0 0 1 1]); 
    if save_to_folder; saveas(gcf, strcat(mydir, meshstr, string(i-1), '_sharp.png')); end
    try
        fprintf(fileID, 'Final Obj: %10.6f', c);
        fclose(fileID);
    catch err
        disp(err.message)
    end
    figure; plot(allenergies, 'o-'); axis tight
    if save_to_folder  
        saveas(gcf, strcat(mydir, meshstr, 'energyplot.png'));
    end
    figure;
    plot(allenergiessame, 'o-'); axis tight
end
