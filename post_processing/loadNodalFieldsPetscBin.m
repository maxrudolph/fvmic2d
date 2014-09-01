function nf = loadNodalFieldsPetscBin(filename)

% filename='loadNodalFields_0_0.petscbin'; % for testing

nf1 = PetscBinaryRead(filename,'cell');
bag = nf1{end};
nf.NX=bag.NX;
nf.NY=bag.NY;
NX=nf.NX;
NY=nf.NY;
nf.elapsedTime=bag.elapsedTime;

%   ierr=  VecView( nodalFields->vx,viewer );CHKERRQ(ierr);         1
%   ierr=  VecView( nodalFields->vy,viewer );CHKERRQ(ierr);         2
%   ierr=  VecView( nodalFields->vz,viewer );CHKERRQ(ierr);         3
%   ierr=  VecView( nodalFields->p,viewer );CHKERRQ(ierr);          4
%   ierr=  VecView( nodalFields->rho,viewer );CHKERRQ(ierr);        5  
%   ierr=  VecView( nodalFields->rhodot,viewer );CHKERRQ(ierr);     6
%   //ierr=  VecView( nodalFields->D,viewer );CHKERRQ(ierr);        x
%   //ierr=  VecView( nodalFields->Ddot,viewer );CHKERRQ(ierr);     x
%   ierr=  VecView( nodalFields->etaN,viewer );CHKERRQ(ierr);       7   
%   ierr=  VecView( nodalFields->etaS,viewer );CHKERRQ(ierr);       8  
%etavx 9 
%etavy 10
%   ierr=  VecView( nodalFields->kThermal,viewer);CHKERRQ(ierr);    11
%   ierr=  VecView( nodalFields->lastT,viewer );CHKERRQ(ierr);      12
%   ierr=  VecView( nodalFields->Cp,viewer);CHKERRQ(ierr);          13
%   ierr=  VecView( nodalFields->muN,viewer);CHKERRQ(ierr);         14
%   ierr=  VecView( nodalFields->muS,viewer);CHKERRQ(ierr);         15

%   muvx 16 
%   muvy 17 
% 
% #ifdef TEXTURE
%   ierr=  VecView( nodalFields->VPTensorB,viewer);CHKERRQ(ierr);   18
%   ierr=  VecView( nodalFields->VPTensorC,viewer);CHKERRQ(ierr);   19
% #endif
%   ierr=  VecView( nodalFields->soxx,viewer);CHKERRQ(ierr);        16
%   ierr=  VecView( nodalFields->soyy,viewer);CHKERRQ(ierr);        17
%   ierr=  VecView( nodalFields->sozz,viewer);CHKERRQ(ierr);        18
%   ierr=  VecView( nodalFields->soxy,viewer);CHKERRQ(ierr);        19
%   ierr=  VecView( nodalFields->soxz,viewer);CHKERRQ(ierr);        20
%   ierr=  VecView( nodalFields->soyz,viewer);CHKERRQ(ierr);        21
%   /* strain rate components */
%   ierr=  VecView( nodalFields->edotxx,viewer);CHKERRQ(ierr);      22
%   ierr=  VecView( nodalFields->edotyy,viewer);CHKERRQ(ierr);      23
%   ierr=  VecView( nodalFields->edotxy,viewer);CHKERRQ(ierr);      24
%   ierr=  VecView( nodalFields->edotxz,viewer);CHKERRQ(ierr);      25
%   ierr=  VecView( nodalFields->edotyz,viewer);CHKERRQ(ierr);      26
% 
%   wxy 27
%   wxz 28
%   wyz 29


%   ierr=  VecView( nodalFields->ha,viewer);CHKERRQ(ierr);          30
% #ifdef TEXTURE
%   ierr = VecView( nodalFields->strainRateResidual,viewer);CHKERRQ(ierr);
% #endif


nf.vx = reshape(nf1{1},[NX NY])';
nf.vy = reshape(nf1{2},[NX NY])';
nf.vz = reshape(nf1{3},[NX NY])';
nf.p = reshape(nf1{4},[NX NY])';
nf.rho = reshape(nf1{5},[NX NY])';
nf.rhodot = reshape(nf1{6},[NX NY])';
nf.etaN = reshape(nf1{7},[NX NY])';
nf.etaS = reshape(nf1{8},[NX NY])';
nf.etavx = reshape(nf1{9},[NX NY])';
nf.etavy = reshape(nf1{10},[NX NY])';
nf.T = reshape(nf1{12},[NX NY])';
if( length(nf1{18}) == NX*NY*4)
    oft=2;
    nf.Nb = reshape(nf1{18},[4 NX NY]);
else
    oft=0;
end

nf.ha = reshape(nf1{23+oft},[NX NY])';
nf.exx = reshape(nf1{20+oft},[NX NY])';
nf.eyy = reshape(nf1{21+oft},[NX NY])';
nf.exy = reshape(nf1{22+oft},[NX NY])';
nf.exz = reshape(nf1{23+oft},[NX NY])';
nf.eyz = reshape(nf1{24+oft},[NX NY])';

nf.sxx = reshape(nf1{16+oft},[NX NY])';
nf.sxy = reshape(nf1{17+oft},[NX NY])';
nf.sxz = reshape(nf1{18+oft},[NX NY])';
nf.syz = reshape(nf1{19+oft},[NX NY])';

nf.wxy = reshape(nf1{25+oft},[NX NY])';
nf.wxz = reshape(nf1{26+oft},[NX NY])';
nf.wyz = reshape(nf1{27+oft},[NX NY])';

nf.muN = reshape(nf1{12},[NX NY])';
nf.muS = reshape(nf1{13},[NX NY])';
if(oft)
nf.srr = reshape(nf1{26+oft},[NX NY])';
end

