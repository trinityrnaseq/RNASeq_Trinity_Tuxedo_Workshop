library(cummeRbund);


cat("\n","Reading cuffdiff output","\n")
cat("cuff = readCufflinks('diff_out');\n")
cuff = readCufflinks('diff_out');

cat("\n","Plotting FPKM density for genes:","\n")
cat("csDensity(genes(cuff))\n");
print(csDensity(genes(cuff)));
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

cat("\n","Plotting scatter for gene fpkms in Sp_log vs. Sp_plat:","\n")
cat("csScatter(genes(cuff), 'Sp_log', 'Sp_plat'))\n");
print(csScatter(genes(cuff), 'Sp_log', 'Sp_plat'));
dev.flush()
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

cat("\n","Plotting scatter matrix:","\n")
cat("csScatterMatrix(genes(cuff)\n");
print(csScatterMatrix(genes(cuff)));
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

cat("\n","Plotting volcano matrix:","\n")
cat("csVolcanoMatrix(genes(cuff), 'Sp_log', 'Sp_plat')\n")
print(csVolcanoMatrix(genes(cuff), 'Sp_log', 'Sp_plat'));
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

cat("\n","Get genes found as diff expressed:","\n")
cat("gene_diff_data = diffData(genes(cuff));\n")
gene_diff_data = diffData(genes(cuff));
# how many genes?
cat("How many genes?\n");
cat("nrow(gene_diff_data)\n");
count = nrow(gene_diff_data);
cat("\n","Have :",count, "\n")
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)


cat("\n","Extract diff expressed genes according to specified criteria:","\n")
cat("sig_gene_data = subset(gene_diff_data,(significant=='yes' | p_value < 0.01))\n");
sig_gene_data = subset(gene_diff_data,(significant=='yes' | p_value < 0.01));
# how many significant diff expr?
cat("How many significant genes selected?\n");
count = nrow(sig_gene_data);
cat("count = nrow(sig_gene_data)\n");
cat("\n", "Found ", count, " genes diff expressed according to criteria: p_value < 0.01", "\n");
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

cat ("\n", "Examples of DE genes:", "\n");
cat("head(sig_gene_data)\n");
print(head(sig_gene_data));
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

ex_gene_id = sig_gene_data$gene_id[1]

cat("Examine the expression values for one of the genes. We'll select" , ex_gene_id, "\n");
cat("ex_gene = getGene(cuff, '", ex_gene_id, "')\n", sep='');
ex_gene = getGene(cuff, ex_gene_id)
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

cat ("Generate the expression plot for this gene.\n");
cat("expressionBarplot(ex_gene, logMode=T, showErrorbars=F)\n");
print(expressionBarplot(ex_gene, logMode=T, showErrorbars=F))
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

cat("Generate a clustered heatmap for all significant genes:\n");
cat("sig_genes = getGenes(cuff, sig_gene_data$gene_id)\n");
sig_genes = getGenes(cuff, sig_gene_data$gene_id);

cat("csHeatmap(sig_genes, cluster='both')\n");
print(csHeatmap(sig_genes, cluster='both'));
cat("\n", "Done with cummeRbund demo", "\n");


