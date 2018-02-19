df = read.table("freq_spectra.txt");
m = as.matrix(df[, -1]);
tmp.m = apply(m, 1, function(x) {x/sum(x)});

#neutral expetation
x = 1:25;
a = 1 / x;
kesai = 1/(x*sum(a));

png("dsi_freq_spectrum.png", width=1200, height=600);
par(cex=1.5);
cols = c("white", "firebrick2", "grey", "black", "yellow", "steelblue2", "tan1");
bp = barplot(t(tmp.m), beside=T, names.arg=c(1:25), col=cols, xlab="Derived allele with size i", ylab="Fraction of SNPs" );
lines(bp[4,], kesai, lty=2);
points(bp[4,], kesai, pch=20);
legend(150,1, legend=c("synonymous", "nonsynonymous", "intron", "intergenic", "dsi-miR-2582b", "dsi-miR-303", "dsi-miR-983"), fill=cols, cex=0.8 );
dev.off();
