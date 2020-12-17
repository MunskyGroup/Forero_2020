load('best_simple_pars_cs.mat')

best_par = parameters;
sub_par = parameters;

sub_par(4) = 5;

[sigred,~,means,SIGc] = get_ac_and_cc_mod_cs(subpar, [0:.2:30] );

[mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap] = load_normalization_variance_gui(1,'G0_intp');

mrna_rnap.sem_cc = mrna_rnap.sem_cc/mrna_rnap.mn_cc(11);
mrna_rnap.mn_cc = mrna_rnap.mn_cc/mrna_rnap.mn_cc(11);

mrna_ser5.sem_cc = mrna_ser5.sem_cc/mrna_ser5.mn_cc(11);
mrna_ser5.mn_cc = mrna_ser5.mn_cc/mrna_ser5.mn_cc(11);

ser5_rnap.sem_cc = ser5_rnap.sem_cc/ser5_rnap.mn_cc(11);
ser5_rnap.mn_cc = ser5_rnap.mn_cc/ser5_rnap.mn_cc(11);

TTc = [-10:10];
TTc = [TTc(1:10) -1:.2:1 TTc(12:end)];


ser5_rnap_sub = [sigred(4,end:-1:2),sigred(2,:)];
s5rn = [ser5_rnap.mn_cc(1:10)'    sigred(4,5:-1:1) 1  sigred(2,1:5) ser5_rnap.mn_cc(12:end)' ];

mrna_ser5_sub = [sigred(7,end:-1:2),sigred(3,:)];
ms5 = [mrna_ser5.mn_cc(1:10)'    sigred(8,5:-1:1) 1  sigred(6,1:5) mrna_ser5.mn_cc(12:end)' ];

mrna_rnap_sub = [sigred(7,end:-1:2),sigred(3,:)];
mrnp = [mrna_rnap.mn_cc(1:10)'    sigred(7,5:-1:1) 1  sigred(3,1:5) mrna_rnap.mn_cc(12:end)' ];

TTa = [1:30];
TTa = [0:.2:1 TTa];

rn = [ rnap.mn_ac(1)' sigred(1,:) rnap.mn_ac(2:end)'  ];
s5 = [ ser5.mn_ac(1)' sigred(5,:) ser5.mn_ac(2:end)'  ];
mrna = [ mrna.mn_ac(1)' sigred(9,:) mrna.mn_ac(2:end)'  ];
