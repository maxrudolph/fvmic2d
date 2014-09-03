function nf = loadNodalFieldsPetscBin2(filename)

% filename='loadNodalFields_0_0.petscbin'; % for testing

nf1 = PetscBinaryRead(filename,'cell',10000);
bag = nf1{end};
nf.NX=bag.NX;
nf.NY=bag.NY;
NX=nf.NX;
NY=nf.NY;
nf.elapsedTime=bag.elapsedTime;

nf.vx = reshape(nf1{1},[NX NY])';
nf.vy = reshape(nf1{2},[NX NY])';
% nf.vz = reshape(nf1{3},[NX NY])';
nf.p = reshape(nf1{3},[NX NY])';
nf.rho = reshape(nf1{4},[NX NY])';
nf.rhodot = reshape(nf1{5},[NX NY])';
nf.etaN = reshape(nf1{6},[NX NY])';
nf.etaS = reshape(nf1{7},[NX NY])';
% nf.etavx = reshape(nf1{9},[NX NY])';
% nf.etavy = reshape(nf1{10},[NX NY])';
nf.T = reshape(nf1{9},[NX NY])';
% if( length(nf1{18}) == NX*NY*4)
%     oft=2;
%     nf.Nb = reshape(nf1{18},[4 NX NY]);
%     nf.Nc = reshape(nf1{19},[4 NX NY]);
% else
%     oft=0;
% end

nf.ha = reshape(nf1{20},[NX NY])';
nf.exx = reshape(nf1{16},[NX NY])';
nf.eyy = reshape(nf1{17},[NX NY])';
nf.exy = reshape(nf1{18},[NX NY])';
% nf.exz = reshape(nf1{27+oft},[NX NY])';
% nf.eyz = reshape(nf1{28+oft},[NX NY])';

nf.sxx = reshape(nf1{13},[NX NY])';
nf.syy = reshape(nf1{14},[NX NY])';
% nf.szz = reshape(nf1{20+oft},[NX NY])';
nf.sxy = reshape(nf1{15},[NX NY])';
% nf.sxz = reshape(nf1{22+oft},[NX NY])';
% nf.syz = reshape(nf1{23+oft},[NX NY])';

nf.wxy = reshape(nf1{19},[NX NY])';
% nf.wxz = reshape(nf1{30+oft},[NX NY])';
% nf.wyz = reshape(nf1{31+oft},[NX NY])';

nf.muN = reshape(nf1{14},[NX NY])';
nf.muS = reshape(nf1{15},[NX NY])';
% nf.muvx = reshape(nf1{16},[NX NY])';
% nf.muvy = reshape(nf1{17},[NX NY])';
% if( length(nf1) >= 34+oft )
%      if(length(nf1{33+oft}) == NX*NY*2)
         tmp = reshape(nf1{21},[2 NX NY]);
         nf.gridx = squeeze(tmp(1,:,:))';
         nf.gridy = squeeze(tmp(2,:,:))';
%      end
% end

% if(oft)
% nf.srr = reshape(nf1{26+oft},[NX NY])';
% end
% nf.nodalHeating = reshape(nf1{end-2},[NX NY])';
