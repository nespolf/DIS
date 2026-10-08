function write_EMC3_imp_input(NAS_EMC3,fnameout,element)
disp('writing input.imp file for EMC3-EIRENE...')
a=readlines('input.imp_31');
b=readlines('input.imp_32_fixed_coeff');
c=readlines('input.imp_33');
d=readlines('input.imp_34');
adas_folder={'/scratch/gpfs/fnespoli/B-atom/'};
writelines(a,fnameout)
if strcmp(element,'B')
    adas_list={'scd89_b.dat','acd89_b.dat','plt89_b.dat','prb89_b.dat'};
elseif strcmp(element,'C')
    adas_list={'scd96_c.dat','acd96_c.dat','plt96_c.dat','prb96_c.dat'};
elseif strcmp(element,'N')
    adas_list={'scd96_n.dat','acd96_n.dat','plt96_n.dat','prb96_n.dat'};
elseif strcmp(element,'Li')
    adas_list={'scd96_li.dat','acd96_li.dat','plt96_li.dat','prb96_li.dat'};
elseif strcmp(element,'Ne')
    adas_list={'scd96_ne.dat','acd96_ne.dat','plt96_ne.dat','prb96_ne.dat'};
else disp('atomic material not supported')
end
for i=1:length(adas_list)
    adas_list{i}=char(append(adas_folder,adas_list{i}));
end
writelines(adas_list,fnameout,WriteMode="append")
writelines(b,fnameout,WriteMode="append")
aa=[size(NAS_EMC3,1) 0];
writelines(num2str(aa),fnameout,WriteMode="append")
writelines(c,fnameout,WriteMode="append")
% for i=1:size(NAS_EMC3,1)
% writelines(num2str(NAS_EMC3(i,:)),fnameout,WriteMode="append")
% end
dlmwrite(fnameout,NAS_EMC3,'delimiter',' ','-append')
writelines(d,fnameout,WriteMode="append")
end