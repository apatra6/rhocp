function pole_figure_mtex_5p11()
    clc; clear;
    
    CS = crystalSymmetry('m-3m', [3.3013 3.3013 3.3013], 'mineral', 'Tantalum');
    % specimen symmetry
    SS = specimenSymmetry('1');
    
    % plotting convention
    setMTEXpref('xAxisDirection','north');
    setMTEXpref('zAxisDirection','outOfPlane');
    fname = 'euler_angle_data.txt';
    
    % plotting
    ori = orientation.load(fname,CS,'ColumnNames',{'phi1','Phi','phi2'});
    odf = calcKernelODF(ori,'halfwidth',10*degree,'resolution',5*degree);
    plotPDF(odf,Miller({0,0,2},{1,1,0},{1,1,1},CS));
    setColorRange([0 2]);

    % saving
    exportgraphics(gcf,'texture.png','Resolution',300);
end