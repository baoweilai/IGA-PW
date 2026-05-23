function  C=Project(Cw)
%Project homogeneous coordinates to physical coordinates.
C=Cw(1:end-1)/Cw(end);
end
