folder1='/p/emc3-eirene/fnespoli/dust_results/temp/';
folder2='/p/emc3-eirene/fnespoli/dust_results/final/';


fname=dir([folder2,'*.mat']);
for i=1:length(fname)
k2=findstr(fname(i).name,'.');
fname(i).name
fname0=fname(i).name(1:(k2(end)-1))
reassamble_dust_results(folder1,folder2,fname0)   
end
    
% fname0='C_n05e+18_rd5_NKU0_theta30';
% reassamble_dust_results(folder1,folder2,fname0)
% fname0='C_n05e+18_rd0.1_NKU0_theta30';
% reassamble_dust_results(folder1,folder2,fname0)
% fname0='C_n05e+18_rd0.8_NKU0_theta30';
% reassamble_dust_results(folder1,folder2,fname0)
