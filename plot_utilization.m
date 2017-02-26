% hw0 = [6329,9829,320,124,2519,4800,1529,5100,111.5,48,79,0.125];
% hw1 = 5536,8178,121,26,2052,3693,1843,4171,53.6,50,106,0.122]
% S.name = {'hw0','hw1'}
clear
S.SliceLUTS = [6329;5536;53200];
S.SliceRegs = [9829;8178;106400];
S.F7Muxes = [320;121;26600]
S.F8Muxes = [124;26;13300]
S.Slice = [2519;2052;13300]
S.LUTasLogic = [4800;3693;53200]
S.LUTasMem = [1529;1843;17400]
S.LUTFFPairs = [5100;4171;53200]
S.BlockRAMTile = [111.5;53.5;140]
S.DSP = [48;50;220]
S.BondedIOB = [79;106;200]
% S.TotPower = [0.125;0.122;0]

T=struct2table(S);
data=table2array(T);
denom = repmat(data(3,:), [2, 1])
barh(data(1:2,:)'./denom'*100);
types = ['hw0','hw1','available'];
names = fieldnames(S);
yticklabels(names)
xlim([0 100])


