function cfg = dg_plot_cfg()
%Return plotting data for DG error checks.
cfg.Example = 'Example_5_3D_DG_PW_IGA_0414';
cfg.fontName = 'Times New Roman';
cfg.figWidth = 4.6;
cfg.figHeight = 3.3;
cfg.margin = struct('left', 0.14, 'right', 0.04, 'bottom', 0.14, 'top', 0.08);
cfg.lineColors = [223 122 094; 060 064 091; 130 178 154; 242 204 142] / 255;
cfg.order1Color = [239 065 067] / 255;
cfg.order15Color = [033 158 188] / 255;
cfg.lineWidth = 2.0;
cfg.orderLineWidth = 1.8;
cfg.markerSize = 7;
end
