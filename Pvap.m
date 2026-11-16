function Pvap=Pvap(Td,Pvap_A,Pvap_B,material)
%Vapor pressure
%Pvap=Pvap(Td,Pvap_A,Pvap_B)
%Pvap=10.^(Pvap_A./Td+Pvap_B); % [Pa] 
Pvap=10.^(Pvap_A./Td+Pvap_B); % [Pa] 

%if strcmp(material,'W')
%	M=xlsread('/home/AI289520/DIS_CODE/DIS/Pvap_W.xlsx');
%	mindiffTd=min(abs(squeeze(M(:,1))-double(Td)));
%	Pvap=M(find(M(:,1)-Td==mindiffTd),2)*1e5;
%end

end

