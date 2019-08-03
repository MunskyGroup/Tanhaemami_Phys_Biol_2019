%% RECORD FILE NAMES IN A TABLE
%Dates at which experiments are performed.
dates = {'072716','081516','072816','081616','080116','081716',...
    '080216','081816','080316','081916','080516','082216',...
    '080616','082316','080716','082616','080816','083016',...
    '081016','090216','081116','091116','081216'}; %Deelta_T's are the same between the two experiment batches??
Alphabet=char('A'+(1:7)-1)';
prefix_un = {''};
for i=1:length(Alphabet)-3
    for j=1:12
        prefix_un{i,j} = [Alphabet(i),num2str(j,'%02i')];
    end
end

prefix_st = {''};
for i=1:length(Alphabet)-4
    for j=1:12
        prefix_st{i,j} = [Alphabet(i+4),num2str(j,'%02i')];
    end
end

rel_dates = [0 19 1 20 5 21 6 22 7 23 9 26 10 27 11 30 12 34 14 37 15 46 16];
Species = {'Pico'};
n=1;
s=1;
for Date = dates %For each of the dates experiments are run.
    for Un_Prefix = prefix_un
        for m=1:size(Un_Prefix,1)
            File_unstained{n,1} = ['../../Data_Files/Data_main/',Un_Prefix{m,1},' ',Date{1},' ',Species{1},' un.fcs.csv']; %Save unstained data in a ".csv" file.
            if exist(File_unstained{n,1},'file')
                FUL{s,1} = File_unstained{n,1};
                s=s+1;
            end
            n=n+1;
        end
    end
end

n=1;
s=1;
for Date = dates %For each of the dates experiments are run.
    for St_Prefix = prefix_st
        for m=1:size(St_Prefix,1)
            File_stained{n,1} = ['../../Data_Files/Data_main/',St_Prefix{m,1},' ',Date{1},' ',Species{1},' st.fcs.csv']; %Save unstained data in a ".csv" file.
            if exist(File_stained{n,1},'file')
                FL{s,1} = File_stained{n,1};
                s=s+1;
            end
            n=n+1;
        end
    end
end