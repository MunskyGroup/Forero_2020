kon = parameters(1);
koff = 1000;
kesc = parameters(3);
kproc = parameters(4);
kin = 1000*parameters(5);
kout = parameters(6);

mu_ctd = kon/(kon+koff)*kin/(kout+kesc)
mu_ctd_on = kin/(kout+kesc)

mu_rna = kon/(kon+koff)*kin/(kout+kesc)*kesc/kproc
mu_rna_on = kin/(kout+kesc)*kesc/kproc

rate_prod = mu_ctd*kesc
rate_prod_on = mu_ctd_on*kesc

burst_refractory_period = 1/kon
burst_duration = 1/koff
pol2_burst_size = kin/koff

burst_efficiency = kesc/(kesc+kout)
mrna_burst_size = pol2_burst_size*burst_efficiency






