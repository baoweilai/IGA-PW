function cfg = paper_style()
%Return fixed plotting style values.
cfg.fontName = 'Times New Roman';
cfg.lineColors = [223 122 094; 060 064 091; 130 178 154; 242 204 142] / 255;
cfg.lineWidth = 2.0;
cfg.markerSize = 7;
cfg.axesFontSize = 10;
cfg.labelSize = 12;
cfg.legendSize = 10;
cfg.axesLineWidth = 1.0;
cfg.tickDir = 'out';
cfg.errorMap = error_colormap(256);
cfg.errorDPI = 600;
end
