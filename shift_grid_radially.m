clear all;
fname1='~/scripts/DIS/plasma_background/LHD/166256_1_pbg_EMC3.mat'
fname2='~/scripts/DIS/plasma_background/LHD/166256_pbg_EMC3_Rax355.mat'
eval(['!cp ',fname1,' ',fname2])
load(fname2,'mesh')
mesh.xcells=mesh.xcells-0.05;
save(fname2,'mesh','-append')