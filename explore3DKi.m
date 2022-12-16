%Code to explore potential ki values for 3D track code with SB estimation
%hard coded to off/perfect information. This will feed into new weight
%caluclations
alg=1;
disp(datetime)
disp('Beginnith the journey into the great unknown')
bayesSBEST=true;

trials=3;
Dall=10e-12;
N=8.75e4;%1.75sec trajs
tau=20e-6;
s=[5,25, 150,300,50,100]*1e3/0.0428;
b=[0,0, 50,200,100,10]*1e3;
cases=[s',b'];
csnum=size(cases);
csnum(2)=[];

% T=table('Size',[csnum, 11],'VariableTypes',["double","double","double","double","double","double","double","double","double","double","double"]);
% T.Properties.VariableNames={'Signal','Background','sbr (s:b)','3D Kalman Ki_xy','3D Kalman Ki_z','2D Bayes + 1D Kalman Ki_xy','2D Bayes + 1D Kalman Ki_z','2D Bayesian + 1D Bayesian Ki_xy','2D Bayesian + 1D Bayesian Ki_z','3D Bayesian Ki_xy','3D Bayesian Ki_z'};
% for i=1:csnum
%     names(i)=strcat("Case ",num2str(i));
% end
% T.Properties.RowNames=names;
% T(:,[1,2])=array2table(cases);
% T(:,3)=array2table(round(cases(:,1)./(cases(:,2))));

%get top level path
pathtop=uigetdir;
cd(pathtop)
datename='221207';

%set title for data save
if alg == 1
    tit = 'Kalman_XY_Z';
    titb = 'Kalman XY, Z';
    chkptaddendum ='3K';
elseif alg == 2
    tit = 'CBWBS_XY_Kal_Z';
    titb = 'COBWEBS XY + Kalman Z';
    chkptaddendum ='2C1K';
elseif alg == 3
    tit = 'CBWBS_XY_Z';
    titb = 'COBWEBS XY, Z';
    chkptaddendum ='2C1C';
else
    error('invalid algorithm choice')
end
oscil=false;
lastchkpt=datetime;

%% run piezo shit
if alg == 1
ki_xy=0:0.025:0.5;
else
    ki_xy=0:0.01:0.3;
end
ki_z=0:0.01:0.1;

resp='galvo_xy_piezo_z';
cd(pathtop)
if bayesSBEST
    subfoldername=[datename '_' tit '_explore_kikp_' resp '_sbest_true'];
else
    subfoldername=[datename '_' tit '_explore_kikp_' resp '_sbest_false'];
end
mkdir(subfoldername);
pathmid=[pathtop '\' subfoldername];
cd(pathmid)


%% Now, see if something bad has happened.
%tell matlab where something bad would be if it did happen or where to save
%temp data.

%If something bad has happened in folder of relevance load progress and flag array.
if isfile([pathmid '/checkpoint_a' chkptaddendum '.MAT'])||isfile([pathmid '/checkpoint_b' chkptaddendum '.MAT'])
    disp('resuming from prior save file')
    FileInfo = dir(['checkpoint_a' chkptaddendum '.MAT']);
    TimeStampA = datetime(FileInfo.date);
    FileInfo = dir(['checkpoint_b' chkptaddendum '.MAT']);
    TimeStampB = datetime(FileInfo.date);
    if TimeStampA>TimeStampB
        load([pathmid '/checkpoint_a' chkptaddendum '.MAT'])
    else
        load([pathmid '/checkpoint_b' chkptaddendum '.MAT'])
    end
else %if nothing bad has happened initialize flag array to run all trials.
    disp('No prior progress exists. Starting from scratch')
    runtrue=true(length(ki_xy),length(ki_z),trials,csnum,length(Dall));
    dur_s=zeros(length(ki_xy),length(ki_z),trials,csnum,length(Dall));
    netmeanerr_um=zeros(length(ki_xy),length(ki_z),trials,csnum,length(Dall));
    xystgdev_m=zeros(length(ki_xy),length(ki_z),trials,csnum,length(Dall));
    zstgdev_m=zeros(length(ki_xy),length(ki_z),trials,csnum,length(Dall));
end

chkptvarnames={'trials','Dall','N','tau','cases','csnum','pathtop','datename','ki_xy','ki_z','pathmid','resp','runtrue','dur_s','netmeanerr_um','xystgdev_m','zstgdev_m'};
outputvarnames={'trials','Dall','N','tau','cases','ki_xy','ki_z','resp','dur_s','netmeanerr_um','xystgdev_m','zstgdev_m','tit','titb'};

if alg ~= 1 % add signal background estimation toggle to saved data
    chkptvarnames=[chkptvarnames,{'bayesSBEST'}];
    outputvarnames=[outputvarnames,{'bayesSBEST'}];
end
%% Begin Run
for q=1:length(Dall)
    disp(['Beginning feedback parameter optimization for ' titb])
    D=Dall(q);
    cd(pathmid)
    
    %generate data
    for r=1:csnum
        s=cases(r,1);
        b=cases(r,2);
        for j=1:length(ki_xy)
            for k=1:length(ki_z)
                for i=1:trials
                    if runtrue(j,k,i,r,q)
                        if alg == 1
                            [trkerr,len,tauout,~,~,stg_out,part_out,~,Nout]=track_Kalman_3D(D,s,b,[],ki_xy(j),ki_z(k),N,tau);
                        elseif alg == 2
                            [trkerr,len,tauout,~,~,stg_out,part_out,~,Nout]=track_XYBayesZKalm(D,s,b,[],ki_xy(j),ki_z(k),N,tau,'sbest',bayesSBEST);
                        elseif alg == 3
                            [trkerr,len,tauout,~,~,stg_out,part_out,~,Nout]=track_XYBayesZBayes(D,s,b,[],ki_xy(j),ki_z(k),N,tau,'sbest',bayesSBEST);
                        end
                        netmeanerr_um(j,k,i,r,q)=mean(sqrt(sum((stg_out-part_out).^2,2)))*1e6;
                        xystgdev_m(j,k,i,r,q)=std(trkerr(1:len,1));
                        zstgdev_m(j,k,i,r,q)=std(trkerr(1:len,2));
                        dur_s(j,k,i,r,q)=double(len)*tauout;
                        runtrue(j,k,i,r,q)=false;
                    end
                    %checkpoint save file
                    if datetime-lastchkpt>duration(2,0,0)%if it has been 2 hours since the last checkpoint, go ahead and do a checkpoint
                        if oscil
                            save(['checkpoint_a' chkptaddendum '.mat'],chkptvarnames{:});
                            oscil=false;
                        else
                            save(['checkpoint_b' chkptaddendum '.mat'],chkptvarnames{:});
                            oscil=true;
                        end
                    end
                end
            end
            disp(['ki xy ' num2str(j) '/' num2str(length(ki_xy)) ' ' datestr(datetime('now'))])
        end
        disp(['CASE ' num2str(r) '/' num2str(csnum) ' FINISHED'])
    end
    disp ('ONE D DONE')
end
save([datename '_' tit '_' resp 'explore_kikp_sb_estimation.MAT'],outputvarnames{:})

delete([pathmid '/checkpoint_a' chkptaddendum '.MAT'])
delete([pathmid '/checkpoint_b' chkptaddendum '.MAT'])

%communicate
disp(datetime)
disp(['Finished feedback parameter exploration for ' titb])