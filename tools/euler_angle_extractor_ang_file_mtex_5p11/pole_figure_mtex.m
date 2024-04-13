function pole_figure_mtex()
%% MATLAB script to plot/reduce input euler angles
    clc;clear;
    
    % crystal symmetry
    CS = crystalSymmetry('m-3m', [4.0478 4.0478 4.0478], 'mineral', 'Aluminium');
    
    % specimen symmetry
    SS = specimenSymmetry('1');
    
    % plotting convention
    setMTEXpref('xAxisDirection','north');
    setMTEXpref('zAxisDirection','outOfPlane');
    
    % file to be imported
    fname = 'full_data_euler.txt';
    
    % plotting
    ori = orientation.load(fname,CS,'ColumnNames',{'phi1','Phi','phi2'});
    odf = calcKernelODF(ori,'halfwidth',10*degree,'resolution',5*degree);
    plotPDF(odf,Miller({0,0,2},{1,1,0},{1,1,1},CS));
    setColorRange([0 2]);

    % saving
    exportgraphics(gcf,'texture.png','Resolution',300);

    % reduce the input texture
    number_of_grains = 512;
    export_VPSC(odf,'reduced_data_euler.txt','points',number_of_grains);
end