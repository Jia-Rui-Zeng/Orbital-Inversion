clear
clc

hopping = -2.7;  hubbard_U = 4; hubbard_V = 0.55;
[ho,loca] = generateHamiltonian(hopping);

num_chain = size(ho,1);
num_chain1 = ceil(num_chain/2)*2;

% guess diagonal terms
h1tm = 0.5*eye(num_chain);
h1tm(1,1) = h1tm(1,1)*2;

% up spin
[wfup,energies_up] = eig(ho+h1tm);
energies_up = diag(energies_up);

% down spin
h2tm = 0.5*eye(num_chain);
[wfdown,energies_down] = eig(ho+h2tm);
energies_down = diag(energies_down);

% begining iteration
err1 = inf;  err2 = inf;
for count = 1:1e2
    h1 = ho;  h2 = ho;
    for jj = 1:num_chain
        % adding on-site repulsion
        h1(jj,jj) = hubbard_U*sum(abs(wfdown(jj,1:num_chain1/2-1)).^2);
        h2(jj,jj) = hubbard_U*sum(abs(wfup(jj,1:num_chain1/2)).^2);
        for kk = 1:num_chain  % adding NN repulsion
            if kk~=jj && ho(jj,kk)~=0 && (jj==10||jj>25)
                h1(jj,jj) = h1(jj,jj) +hubbard_V*sum(abs(wfup(kk,1:num_chain1/2)).^2);
                h2(jj,jj) = h2(jj,jj) +hubbard_V*sum(abs(wfdown(kk,1:num_chain1/2-1)).^2);
            end
        end
    end
    [wfup,energies_up1] = eig(h1);
    [wfdown,energies_down1] = eig(h2);
    energies_up1 = diag(energies_up1);
    energies_down1 = diag(energies_down1);
    err1 = sqrt(mean((energies_up1-energies_up).^2));
    err2 = sqrt(mean((energies_down1-energies_down).^2));
    energies_up = energies_up1;
    energies_down = energies_down1;
end
energies = [energies_up, energies_down];
clear h1tm h2tm
clear energies_up1 energies_down1 energies_up energies_down

%% drawing molecular orbitals
for ii = num_chain1/2-2-1:num_chain1/2+2
    figure
    vec = wfup(:,ii);
    % drawing carbon skeletons
    for jj = 1:length(vec)
        for kk = jj+1:length(vec)
            length_tm = sqrt(sum((loca(jj,:)-loca(kk,:)).^2));
            if length_tm>0.9 && length_tm<1.1
                plot(loca([jj,kk],1)+2,loca([jj,kk],2)+2.5,'k','LineWidth',1)
                hold on
                plot(loca([jj,kk],1)-2,loca([jj,kk],2)-2.5,'k','LineWidth',1)
            end
        end
    end
    % drawing up-spin orbitals
    for jj = 1:length(vec)
        if vec(jj) >1e-6
            scatter(loca(jj,1)-2,loca(jj,2)-2.5,vec(jj)*500, ...
                'MarkerEdgeColor',[0,176,176]/256,'LineWidth',1)
        elseif vec(jj) <-1e-6
            scatter(loca(jj,1)-2,loca(jj,2)-2.5,-vec(jj)*500,'filled', ...
                'MarkerEdgeColor',[0,176,176]/256,'MarkerFaceColor',[0,176,176]/256)
        end
    end
    % drawing down-spin orbitals
    vec = wfdown(:,ii);
    for jj = 1:length(vec)
        if vec(jj) >1e-6
            scatter(loca(jj,1)+2,loca(jj,2)+2.5,vec(jj)*500, ...
                'MarkerEdgeColor',[0,176,176]/256,'LineWidth',1)
        elseif vec(jj) <-1e-6
            scatter(loca(jj,1)+2,loca(jj,2)+2.5,-vec(jj)*500,'filled', ...
                'MarkerEdgeColor',[0,176,176]/256,'MarkerFaceColor',[0,176,176]/256)
        end
    end
end

%%
function [hamiltonian,loca] = generateHamiltonian(hopping)
% 
num_benzen = 2;
ho = zeros(6);
for jj = 1:5
    ho(jj,jj+1) = hopping;
end
ho(1,6) = hopping;
h1 = ho;
if num_benzen >1
    for jj = 1:num_benzen-1
        h1 = blkdiag(h1,ho);
        h1(end-6-2,end-5) = hopping*cos(34.8/180*pi);
    end
end
loca1 = [0,0; 0.5,-0.86603; 1.5,-0.86603;
    2,0; 1.5,0.86603; 0.5,0.86603];
loca = loca1;
if num_benzen >1
    for jj = 1:num_benzen-1
        locatm = loca1;
        locatm = locatm +loca(end-2,:) +[1,0];
        loca = [loca; locatm];
    end
end

h2 = zeros(13);
for ii = 1:12
    h2(ii,ii+1) = hopping;
end
h2(1,6) = hopping;
h2(8,13) = hopping;
h2(6,7) = hopping*cos(48.3/180*pi);  h2(7,8) = hopping*cos(48.3/180*pi);
loca2 = [-1.5,-0.86603; -2,-1.73206; -1.5,-2.59809; -0.5,-2.59809; 0,-1.73206; -0.5,-0.86603;
    0,0;
    -0.5,0.86603; 0,1.73206; -0.5,2.59809; -1.5,2.59809; -2,1.73206; -1.5,0.86603];

hamiltonian = blkdiag(h1,h2);
hamiltonian(1,size(h1,1)+7) = hopping;
hamiltonian(1,size(h1,1)+7) = hopping*cos(47.0/180*pi);

htm = zeros(6);
htm(1,2) = hopping; htm(5,6) = hopping;
htm(1,3) = hopping; htm(3,5) = hopping;
htm(2,4) = hopping; htm(4,6) = hopping;
h1tm = htm;
for jj = 1:2
    h1tm = blkdiag(h1tm,zeros(4));
    h1tm(end-5:end,end-5:end) = htm;
end
loca1tm = [0,-0.5; 0,0.5; 0.86603,-1;
    0.86603,1; 1.73206,-0.5; 1.73206,0.5];
loca3tm = loca1tm;
for jj = 1:2
    loca2tm = loca1tm;
    loca2tm = loca2tm +[jj*1.73206,0];
    loca3tm = [loca3tm; loca2tm];
    if num_benzen >1; loca3tm(end-5:end-4,:) = []; end
end
loca4tm = loca3tm;
loca3tm(:,1) = loca4tm(:,2);
loca3tm(:,2) = loca4tm(:,1);
hamiltonian = blkdiag(hamiltonian,h1tm);
hamiltonian(size(h1,1)-2,size(h1,1)+13+size(h1tm)/2) = hopping*cos(70.2/180*pi);

hamiltonian = hamiltonian +hamiltonian';

loca = loca +loca2(7,:) +[1,0];
locatm = loca2;
locatm(:,1) = -locatm(:,1);
locatm = locatm -locatm(7,:) +loca(end-2,:) +[1,0];
loca = [loca; loca2];
loca3tm = loca3tm -loca3tm(7,:) +loca(size(h1,1)-2,:) +[1,0]; loca = [loca; loca3tm];
loca = loca -[mean(loca(:,1)),0];

end

