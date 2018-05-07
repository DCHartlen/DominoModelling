function Data = glstat_parser(filename)

% Parse GLSTAT file from LS-DYNA output into MATLAB
% Returns table of GLSTAT data
%
% Created By:     D. Hartlen
% Date:           10-May-2016
% Updated By:
% Date:
%


% Initialization

% Parameters
n_block = 0;
Data = struct([]);

% Open file
fid = fopen(filename);

% skip header information
for i=1:4
    fgetl(fid);
end

% Begin reading in data
while(~feof(fid))
    % Read in the data
    n_block = n_block+1;
    block_1 = textscan(fid,'%34s%12f',18,'WhiteSpace','.','HeaderLines',1);
    block_2 = textscan(fid,'%34s(nanosec)%12f',1,'WhiteSpace','.');
    fgetl(fid);fgetl(fid);fgetl(fid); %skip three lines to next datablock
    % Sort the data
    Data(n_block).('Time') = block_1{1,2}(1);
    Data(n_block).('TimeStep') = block_1{1,2}(2);
    Data(n_block).('KineticEng') = block_1{1,2}(3);
    Data(n_block).('InternalEng') = block_1{1,2}(4);
    Data(n_block).('SpringEng') = block_1{1,2}(5);
%     Data(n_block).('HourglassEng') = block_1{1,2}(6);
    Data(n_block).('DampingEng') = block_1{1,2}(6);
    Data(n_block).('InterfaceEng') = block_1{1,2}(7);
    Data(n_block).('ExternalWork') = block_1{1,2}(8);
    Data(n_block).('ErdKineticEng') = block_1{1,2}(9);
    Data(n_block).('ErdInternalEng') = block_1{1,2}(10);
    Data(n_block).('ErdHourglassEng') = block_1{1,2}(11);
    Data(n_block).('TotalEng') = block_1{1,2}(12);
    Data(n_block).('TotalvInitial') = block_1{1,2}(13);
    Data(n_block).('TotalvInitialNoErd') = block_1{1,2}(14);
    Data(n_block).('XVel') = block_1{1,2}(15);
    Data(n_block).('YVel') = block_1{1,2}(16);
    Data(n_block).('ZVel') = block_1{1,2}(17);
%     Data(n_block).('TimePerZoneCycle') = block_2{1,2}(1);
end

% convert structure to table
Data = struct2table(Data);
fclose(fid);
end





