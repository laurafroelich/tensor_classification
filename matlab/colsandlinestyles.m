
allmethodnames = {'CMDA', 'DATER',...
    'ManTDA', 'ManPDA', ...
    'ManTDA\_sr', 'ManPDA\_sr',...
    'Tucker', 'Parafac',...
    'DGTDA', 'Man\_alternatingmodes', 'Parafac2', ...
    'CMDAManTDA', 'CMDAManPDA', 'CMDAManTDA\_sr', 'CMDAManPDA\_sr', ...
    'BDCA\_dampnewton', 'BDCA', 'BDCATucker', 'Tucker2', 'DATEReig'};
colororder = [
    0.00  0.00  1.00
    0.00  0.50  0.00
    1.00  0.00  0.00
    0.00  0.75  0.75
    0.75  0.00  0.75
    0.75  0.75  0.00
    0.25  0.25  0.25
    0.75  0.25  0.25
    0.95  0.95  0.00
    0.25  0.25  0.75
    0.75  0.75  0.75
    0.00  1.00  0.00
    0.76  0.57  0.17
    0.54  0.63  0.22
    0.34  0.57  0.92
    1.00  0.10  0.60
    0.88  0.75  0.73
    0.10  0.49  0.47
    0.66  0.34  0.65
    0.99  0.41  0.23
    ];
CMDAcol = colororder(1,:);
DATERcol = colororder(2,:);
ManTDAcol = colororder(3,:);
ManPDAcol = colororder(4,:);
ManTDA_normsratiocol = colororder(15,:);
ManPDA_normsratiocol = [0 1 0];
Tuckercol = colororder(14,:);
Parafaccol = colororder(10,:);
DGTDAcol = colororder(5,:);
ManAltcol = colororder(12,:);
Parafac2col = colororder(13,:);
BDCAcol = colororder(14,:);
BDCA_BFGScol = colororder(6,:);
BDCAtucker_BFGScol = colororder(7,:);
shadefact = 0.2;



% http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
colororder = [
    0.00  0.00  0.00
    0.00  73  73
    0  146  146
    255  109  182
    255  182  119
    73  0  146
    0  109  219
    182  109  255
    109  182  255
    182  219  255
    146  0  0
    146  73  0
    219  209  0
    36  255  36
    255  255  109
    ];
colororder = colororder/255;%./repmat(sum(colororder,2), 1, 3);
CMDAcol = colororder(2,:);
DATERcol = colororder(3,:);
ManTDAcol = colororder(4,:);
ManPDAcol = colororder(5,:);
ManTDA_normsratiocol = colororder(6,:);
ManPDA_normsratiocol = colororder(7,:);
Tuckercol = colororder(8,:);
Parafaccol = colororder(9,:);
DGTDAcol = [0 1 0.75];
ManAltcol = colororder(11,:);
Parafac2col = colororder(12,:);
BDCAcol = colororder(13,:);
BDCA_BFGScol = colororder(15,:);
BDCAtucker_BFGScol = colororder(14,:);
Tucker2col = [0 0.8 0.2];
DATEReigcol = [1 0 0];
hfact = 1;
sfact = 0.5;
vfact = 1;
CMDAManTDAcol = hsv2rgb(rgb2hsv(ManTDAcol).*[hfact sfact vfact]);
CMDAManPDAcol = hsv2rgb(rgb2hsv(ManPDAcol).*[hfact sfact vfact]);
CMDAManTDA_normsratiocol = hsv2rgb(rgb2hsv(ManTDA_normsratiocol).*[hfact sfact vfact]);
CMDAManPDA_normsratiocol = hsv2rgb(rgb2hsv(ManPDA_normsratiocol).*[hfact sfact vfact]);

cols = [CMDAcol; DATERcol; ManTDAcol; ManPDAcol;...
    ManTDA_normsratiocol; ManPDA_normsratiocol; Tuckercol; ...
    Parafaccol; DGTDAcol; ManAltcol; Parafac2col; CMDAManTDAcol; CMDAManPDAcol;...
    CMDAManTDA_normsratiocol;...
    CMDAManPDA_normsratiocol;...
    BDCAcol; BDCA_BFGScol; BDCAtucker_BFGScol; Tucker2col; DATEReigcol];

linestyles = [repmat({'-'},1,6), repmat({':'},1,2), '-', 'a', ':',...
    repmat({'-'},1,4), repmat({'--'},1,3), ':', '-'];

gcafontsize = 26;
labelfontsize = 30;
titlefontsize = 36;
legfontsize = 22;