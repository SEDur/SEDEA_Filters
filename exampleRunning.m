%Generate some filter coefs and do some plotting

figure();
filc = Sedea_Rbj_Matlabfilters(5000, 48000, 2.0, -2.0);
somenums = sedea_rbjM_lpf(filc);
subplot(2,2,1);
zplane(somenums(1,:),somenums(2,:));
title('lpf coefs');
legend();
somenums1 = sedea_rbjM_hpf(filc);
subplot(2,2,2);
zplane(somenums1(1,:),somenums1(2,:));
title('hpf coefs');
legend();
somenums2 = sedea_rbjM_bpfcq(filc);
subplot(2,2,3);
zplane(somenums2(1,:),somenums1(2,:));
title('bpf coefs');
legend();
somenums3 = sedea_rbjM_bpfcg(filc);
somenums4 = sedea_rbjM_notch(filc);
subplot(2,2,4);
zplane(somenums4(1,:),somenums4(2,:));
title('notch coefs');
legend();
somenums5 = sedea_rbjM_apf(filc);
somenums6 = sedea_rbjM_pek(filc);
somenums7 = sedea_rbjM_ls(filc);
somenums8 = sedea_rbjM_hs(filc);

%generate a noise signal

%Filter the noise

%plot the outputs
figure();
