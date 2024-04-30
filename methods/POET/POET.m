function [Sy,Su,F,B]=POET(Return,funpath,Rpath)

sep=filesep;

RscriptFileName=[funpath,sep,'POET_R.R'];
save([funpath,sep,'data_temp.mat'],'Return')

commandline=['"',Rpath,sep,'R.exe" CMD BATCH "',RscriptFileName,'"'];
system(commandline);

load([funpath,sep,'data_temp.mat'],'Sy','Su','F','B')
delete([funpath,sep,'data_temp.mat'])
delete([funpath,sep,'POET_R.Rout'])
delete([funpath,sep,'.Rdata'])
end