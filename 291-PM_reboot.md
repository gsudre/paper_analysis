# 2022-05-17 11:29:00

The goal for this note is to summarize all the code used in the PM paper into
one page, so that it can easily be shared later. Also, I'll run some of the
other analysis suggested by the reviewers... almost there!

# Main analysis

## DGE

```r
pca_DGE_clean = function(myregion, fm_str, varvec, wnhOnly=FALSE, loo=0) {
    # varvec is a list of variable name and the groups to keep
    data = read.table('~/data/post_mortem_adhd/adhd_rnaseq_counts.txt', header=1)
    rownames(data) = data[,1]
    data[,1] = NULL
    data = round(data)
    sub_name = gsub(x=colnames(data), pattern='X', replacement='')
    colnames(data) = sub_name
    # this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
    # its other replicate
    data = data[, ! colnames(data) %in% c('66552')]
    # outliers based on PCA plots
    outliers = c('68080','68096', '68108', '68084', '68082')
    data = data[, ! colnames(data) %in% outliers]

    library(gdata)
    df = read.xls('~/data/post_mortem_adhd/POST_MORTEM_META_DATA_JAN_2021.xlsx')
    if (!is.na(varvec)) {
        keep_me = c()
        for (key in 1:length(varvec)) {
            for (val in 1:length(varvec[[key]])) {
                keep_me = c(keep_me,
                            which(df[, names(varvec)[key]] == varvec[[key]][val]))
            }
        }
        df = df[unique(keep_me), ]
    }

    # reduce sample to White non-Hispanics only if requested
    if (wnhOnly) {
        imWNH = which(df$C1 > 0 & df$C2 < -.075)
        df = df[imWNH, ]
    }

    # sync data and df
    data = data[, colnames(data) %in% df$submitted_name]
    df = df[df$submitted_name %in% colnames(data), ]
    df = df[order(df$submitted_name), ]
    data = data[, order(df$submitted_name)]

    # keep only data for either ACC or Caudate
    keep_me = df$Region == myregion
    data = data[, keep_me]
    df = df[keep_me, ]

    # leave out subject if needed
    if (loo > 0) {
        cat('leaving out subject', loo, '\n')
        data = data[, -loo]
        df = df[-loo,]
        # rename samples so DESeq2 doesn't complain later
        cnames = sapply(1:ncol(data), function(x) sprintf('sample%d', x))
        colnames(data) = cnames
        rownames(df) = cnames
    }

    # cleaning up some variables
    df$Individual = factor(df$hbcc_brain_id)
    df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
    df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
    df$MoD = factor(df$Manner.of.Death)
    df$Sex = factor(df$Sex)
    df$batch = factor(df$batch)
    df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
    df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
    df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
    df$SUB2 = 'no'
    df[df$substance_group > 0, 'SUB2'] = 'yes'
    df$SUB2 = factor(df$SUB2)
    df$substance_group = factor(df$substance_group)
    df$comorbid_group = factor(df$comorbid_group_update)
    df$evidence_level = factor(df$evidence_level)
    df$brainbank = factor(df$bainbank)
    # replace the one subject missing population PCs by the median of their
    # self-declared race and ethnicity
    idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
    pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
    med_pop = apply(df[idx, pop_pcs], 2, median)
    df[which(is.na(df$C1)), pop_pcs] = med_pop
    # combine brain bank and batch variables
    df$BBB = factor(sapply(1:nrow(df),
                            function(x) sprintf('%s_%s',
                                        as.character(df[x,'brainbank']),
                                        as.character(df[x, 'batch']))))
    df$BBB2 = NA                                                                        
    df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
    df[df$batch==3, 'BBB2'] = 2                                                         
    df[df$batch==4, 'BBB2'] = 3      
    df$BBB2 = factor(df$BBB2)
    imWNH = which(df$C1 > 0 & df$C2 < -.075)
    df$POP_BIN = 'other'
    df[imWNH, 'POP_BIN'] = 'WNH'
    df$POP_BIN = factor(df$POP_BIN)  

    # removing non-autosome genes
    library(GenomicFeatures)
    txdb <- loadDb('~/data/post_mortem_adhd/Homo_sapies.GRCh38.97.sqlite')
    txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
                "GENEID")
    bt = read.csv('~/data/post_mortem_adhd/Homo_sapiens.GRCh38.97_biotypes.csv')
    bt_slim = bt[, c('gene_id', 'gene_biotype')]
    bt_slim = bt_slim[!duplicated(bt_slim),]
    txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
    tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
    tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
    imautosome = which(tx_meta$TXCHROM != 'X' &
                    tx_meta$TXCHROM != 'Y' &
                    tx_meta$TXCHROM != 'MT')
    data = data[imautosome, ]
    tx_meta = tx_meta[imautosome, ]

    # remove constant genes (including zeros) as it breaks PCA
    const_genes = apply(data, 1, sd) == 0
    data = data[!const_genes, ]

    library("DESeq2")
    # making sure any numeric covariates are scaled
    num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
            'C1', 'C2', 'C3', 'C4', 'C5')
    for (var in num_vars) {
        df[, var] = scale(df[, var])
    }

    # creating DESeq2 object
    cat('Running', fm_str, '\n')
    dds <- DESeqDataSetFromMatrix(countData = data,
                                  colData = df,
                                  design = as.formula(fm_str))
    # remove genes based on how many subjects have zero counts
    # the minimum number of subjects with zero counts in a gene for the gene
    # to be kept corresponds to the smallest group size (number of Cases)
    min_subjs = min(table(df$Diagnosis))
    keep <- rowSums(counts(dds) == 0) <= min_subjs
    dds <- dds[keep,]

    # remove genes based on filterByExpr()
    library(edgeR)
    design = model.matrix(as.formula(fm_str), data=colData(dds))
    isexpr <- filterByExpr(counts(dds), design=design)
    ddsExpr = dds[isexpr, ]

    # run the actual DESeq2 analysis
    ddsExpr = DESeq(ddsExpr)

    return(ddsExpr)
}

g = 'main'
varvec = NA
fm_str = '~ RINe + C1 + BBB2 + comorbid_group + SUB2 + Diagnosis'
dds.ACC = pca_DGE_clean('ACC', fm_str, varvec)
dds.Caudate = pca_DGE_clean('Caudate', fm_str, varvec)
save(dds.ACC, dds.Caudate,
     file=sprintf('~/data/post_mortem_adhd/results/%s_DGE.RData', g))
```

Export main DGE results to CSV:

```r
library(DESeq2)
mydir = '~/data/post_mortem_adhd/'
mart = readRDS(sprintf('%s/mart_rnaseq.rds', mydir))

library(GenomicFeatures)
txdb <- loadDb(sprintf('%s/Homo_sapies.GRCh38.97.sqlite', mydir))
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
bt = read.csv(sprintf('%s/Homo_sapiens.GRCh38.97_biotypes.csv', mydir))
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]

load(sprintf('%s/results/%s_DGE.RData', mydir, g))
for (r in c('ACC', 'Caudate')) {
    res_str = sprintf('res = results(dds.%s, name = "Diagnosis_Case_vs_Control", alpha=.05)',
                    r)
    eval(parse(text=res_str))
    fname = sprintf('%s/results/DGE_%s_%s_annot.csv', mydir, r, g)

    df = as.data.frame(res)
    colnames(df)[ncol(df)] = 'padj.FDR'
    # grabbing HUGO IDs
    df$GENEID = substr(rownames(df), 1, 15)
    df2 = merge(df, mart, sort=F,
                by.x='GENEID', by.y='ensembl_gene_id', all.x=T, all.y=F)
    df2 = merge(df2, bt_slim, sort=F,
                by.x='GENEID', by.y='gene_id', all.x=T, all.y=F)
    df2 = df2[order(df2$pvalue), ]
    
    write.csv(df2, row.names=F, file=fname)
}
```

## GSEA

Start with developmental sets. They don't have built-in DBs, so we need to use
our own file. First, we create our GMT:

```r
library(ActivePathways)
library(ABAEnrichment)

dev_names = c('prenatal', 'infant (0-2 yrs)', 'child (3-11 yrs)',
              'adolescent (12-19 yrs)', 'adult (>19 yrs)')
# just getting a template for GMT
gmt = read.GMT('~/data/post_mortem_adhd/hsapiens_disease_Disgenet_entrezgene.gmt')

co = .90

for (r in c('ACC', 'Caudate')) {
    struc_id = ifelse(r == 'ACC', 'Allen:10278', 'Allen:10333')
    other_id = ifelse(r == 'Caudate', 'Allen:10278', 'Allen:10333')
    anno = get_annotated_genes(structure_ids=struc_id, dataset='5_stages', 
                                cutoff_quantiles=c(co))
    anno2 = get_annotated_genes(structure_ids=other_id, dataset='5_stages', 
                                cutoff_quantiles=c(co))
    other_genes = unique(anno2[, 'anno_gene'])

    junk = gmt[1:6]
    # calculating genes present in all dev stages
    idx = anno$age_category==1 & anno$cutoff==co
    genes_overlap = unique(anno[idx, 'anno_gene'])
    for (s in 2:5) {
        idx = anno$age_category==s & anno$cutoff==co
        g2 = unique(anno[idx, 'anno_gene'])
        genes_overlap = intersect(genes_overlap, g2)
    }
    # calculating genes unique to specific dev stages
    genes_unique = list()
    for (s in 1:5) {
        others = setdiff(1:5, s)
        idx = anno$age_category==s & anno$cutoff==co
        g = unique(anno[idx, 'anno_gene'])
        for (s2 in others) {
            idx = anno$age_category==s2 & anno$cutoff==co
            g2 = unique(anno[idx, 'anno_gene'])
            rm_me = g %in% g2
            g = g[!rm_me]
        }
        genes_unique[[s]] = unique(g)
    }
    g = genes_overlap
    # first set is overlap genes only present in this current region
    cnt = 1
    g = g[! g %in% other_genes]
    tmp = list(id = 'pandev', genes = unique(g), name = 'pandev')
    junk[[cnt]] = tmp
    cnt = cnt + 1
    # in each dev stage
    for (i in 1:5) {
        # we only use genes unique to that dev stage
        g = genes_unique[[i]]
        # that are not expressed in the other brain region
        g = g[! g %in% other_genes]
        tmp = list(id = sprintf('dev%d', i),
                    genes = unique(g), name = dev_names[i])
        junk[[cnt]] = tmp
        cnt = cnt + 1
    }
    gmt_name = sprintf('~/data/post_mortem_adhd/%s_dev_sets.gmt', r)
    write.GMT(junk, gmt_name)
}
```

Now we can run GSEA:

```r
library(WebGestaltR)
library(DESeq2)

data_dir = '~/data/post_mortem_adhd/results/'
ncpu=31

g = 'main'

load(sprintf('%s/%s_DGE.RData', data_dir, g))
for (region in c('ACC', 'Caudate')) {
    res_str = sprintf('dds = dds.%s', region)
    eval(parse(text=res_str))

    res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))
    
    ranks = -log10(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)
    
    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    res_str = sprintf('GSEA_%s_%s_developmental', g, region)

    db_file = sprintf('%s/../%s_dev_sets.gmt', data_dir, region)
    set.seed(42)
    cat(db_file, '\n')
    project_name = sprintf('%s_10K', res_str)
    enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                        organism="hsapiens",
                        enrichDatabaseFile=db_file,
                        enrichDatabaseType="genesymbol",
                        interestGene=tmp2,
                        outputDirectory = data_dir,
                        interestGeneType="ensembl_gene_id",
                        sigMethod="top", topThr=50,
                        minNum=3, projectName=project_name,
                        isOutput=F, isParallel=T,
                        nThreads=ncpu, perNum=10000, maxNum=2000))
    out_fname = sprintf('%s/%s.csv', data_dir, project_name)
    write.csv(enrichResult, file=out_fname, row.names=F, quote=T)
}
```

And then we do the gene ontology sets, just using the WebGestalt database:

```r
library(WebGestaltR)
library(DESeq2)

data_dir = '~/data/post_mortem_adhd/results'
ncpu = 31
region = 'ACC'  # ACC or Caudate

g = 'main'
load(sprintf('%s/%s_DGE.RData', data_dir, g))

res_str = sprintf('dds = dds.%s', region)
eval(parse(text=res_str))
res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))

ranks = -log10(res$pvalue) * sign(res$log2FoldChange)
geneid = substring(rownames(res), 1, 15)

tmp2 = data.frame(geneid=geneid, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

res_str = sprintf('GSEA_%s_%s', g, region)
DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    set.seed(42)
    cat(res_str, db, '\n')
    project_name = sprintf('%s_%s_10K', res_str, db)
    enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="ensembl_gene_id",
                                sigMethod="top", topThr=150000,
                                outputDirectory = data_dir,
                                minNum=3, projectName=project_name,
                                isOutput=F, isParallel=T,
                                nThreads=ncpu, perNum=10000))
    out_fname = sprintf('%s/%s.csv', data_dir, project_name)
    write.csv(enrichResult, file=out_fname, row.names=F, quote=T)
}
```

And finally we do GSEA for CNV-implicated genes. Like the develoepmental
analysis, we need to create the GMT file for WebGestalt:

```r
library(ActivePathways)
# just getting a template
gmt = read.GMT('~/data/post_mortem_adhd/ACC_dev_sets.gmt')
# also available in supplemental table
mygenes = c('ASTN2', 'TRIM32', 'BCHE', 'CHRNA7', 'COPA', 'CPSF2', 'CSMD1',
            'DMRTB1', 'GRM1', 'GRM5', 'GRM7', 'GRM8', 'MAPK1', 'MSH3', 'MYC',
            'NCSTN', 'NDUFA5', 'NDUFB1', 'NRXN1', 'PARK2', 'PDCD6IP', 'PEA15',
            'PLOD3', 'POLR1A', 'POLR2B', 'POLR3C', 'PPM1F', 'RBFOX1', 'SEPT5',
            'SLC2A3', 'TBL1XR1', 'TTC27', 'TUBA3C', 'TUBGCP2', 'WASL', 'WWOX')
tmp = list(id = 'CNV_conservative', genes = mygenes, name = 'CNV')
gmt[[1]] = tmp
mygenes = c('COPA', 'CPSF2', 'CSMD1', 'DMRTB1', 'MAPK1', 'MSH3', 'MYC', 'NCSTN',
            'NDUFA5', 'NDUFB1', 'PARK2', 'PDCD6IP', 'PEA15', 'PLOD3', 'POLR1A',
            'POLR2B', 'POLR3C', 'PPM1F', 'RBFOX1', 'SEPT5', 'TBL1XR1', 'TTC27',
            'TUBA3C', 'TUBGCP2', 'WASL', 'WWOX')
tmp = list(id = 'CNV_Harich', genes = mygenes, name = 'CNV2')
gmt[[2]] = tmp
# doing more than one set so WebGestalt won't complain
gmt[3:length(gmt)] = NULL
gmt_name = '~/data/post_mortem_adhd/CNV_gene_sets.gmt'
write.GMT(gmt, gmt_name)
```

Now we run it:

```r
library(WebGestaltR)
library(DESeq2)

data_dir = '~/data/post_mortem_adhd/results/'
ncpu=24

g = 'main'

load(sprintf('%s/%s_DGE.RData', data_dir, g))
for (region in c('ACC', 'Caudate')) {
    res_str = sprintf('dds = dds.%s', region)
    eval(parse(text=res_str))

    res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))
    
    ranks = -log(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)
    
    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    res_str = sprintf('GSEA_%s_%s', g, region)

    DBs = c('CNV_gene_sets')
    for (db in DBs) {
        set.seed(42)
        cat(res_str, db, '\n')
        db_file = sprintf('~/data/post_mortem_adhd/%s.gmt', db)
        project_name = sprintf('%s_%s_10K', res_str, db)
        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                            organism="hsapiens",
                            enrichDatabaseFile=db_file,
                            enrichDatabaseType="genesymbol",
                            interestGene=tmp2,
                            outputDirectory = data_dir,
                            interestGeneType="ensembl_gene_id",
                            sigMethod="top", topThr=50,
                            minNum=3, projectName=project_name,
                            isOutput=F, isParallel=T,
                            nThreads=ncpu, perNum=10000, maxNum=2000))
        out_fname = sprintf('%s/%s.csv', data_dir, project_name)
        write.csv(enrichResult, file=out_fname, row.names=F, quote=T)
    }
}
```


## MAGMA

I downloaded the EUR and AFR 1KG reference files from their website
(https://ctg.cncr.nl/software/magma). Then, we run MAGMA in the command line:

### Generic files

```bash
cd ~/data/post_mortem_adhd/results;
mkdir MAGMA;
cd MAGMA;
module load plink/1.9.0-beta4.4;
echo "g1000_afr.bed g1000_afr.bim g1000_afr.fam" > merge_list.txt;
# they suggest using MAF when subsampling, so let's use it for concatenating too
plink --bfile g1000_eur --merge-list merge_list.txt --maf 1e-5 \
    --flip-scan --make-bed --out g1000_BW;
plink --bfile g1000_eur --exclude g1000_BW-merge.missnp \
    --make-bed --out d1;
plink --bfile g1000_afr --exclude g1000_BW-merge.missnp \
    --make-bed --out d2;
echo "d2.bed d2.bim d2.fam" > merge_list.txt;
plink --bfile d1 --merge-list merge_list.txt --maf 1e-5 \
    --make-bed --out g1000_BW;

module load MAGMA/1.09a;

# annotation step
magma --annotate --seed 42 --snp-loc g1000_BW.bim \
    --gene-loc /usr/local/apps/MAGMA/gene_location/NCBI37.3/NCBI37.3.gene.loc \
    --out annot_BW;
```

Now that we have the generic files MAGMA will need, let's create the gene
covariance file coming from our DEG analysis:

```r
library(org.Hs.eg.db)
library(GenomicFeatures)
G_list0 = select(org.Hs.eg.db, keys(org.Hs.eg.db, 'ENSEMBL'),
                 columns=c('ENSEMBL', 'ENTREZID', 'SYMBOL'), 'ENSEMBL')
library(dplyr)
library(DESeq2)

g = 'main'

load(sprintf('~/data/post_mortem_adhd/results/%s_DGE.RData', g))
for (r in c('ACC', 'Caudate')) {
    res_str = sprintf('dds = dds.%s', r)
    eval(parse(text=res_str))

    res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))

    res$GENEID = substr(rownames(res), 1, 15)
    G_list <- G_list0[!is.na(G_list0$ENSEMBL),]
    G_list = G_list[G_list$ENSEMBL!='',]
    G_list <- G_list[!duplicated(G_list$ENSEMBL),]
    imnamed = res$GENEID %in% G_list$ENSEMBL
    res = res[imnamed, ]
    res2 = merge(res, G_list, sort=F, all.x=F, all.y=F, by.x='GENEID',
                by.y='ENSEMBL')
    ranks = res2 %>% group_by(ENTREZID) %>% slice_min(n=1, pvalue, with_ties=F)
    myres = data.frame(gene=ranks$ENTREZID,
                       rank=sign(ranks$log2FoldChange)*-log10(ranks$pvalue))
    out_fname = sprintf('~/data/post_mortem_adhd/results/MAGMA_%s_dge_%s.tab',
                        g, r)
    write.table(myres, row.names=F, sep='\t', file=out_fname, quote=F)
}
```

### ADHD

It has EUR and multi-ethnic GWAS results. For the main analysis we only use the latter.

```bash
g=main;
magma --bfile g1000_BW --seed 42 --pval ../../adhd_jul2017 N=55374 \
    --gene-annot annot_BW.genes.annot --out genes_ADHD_BW;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_ADHD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_ADHD_${r};
done;
```

### ASD

I only found GWAS for WNH population only. So I'll have to use the _eur BIM file
for the main analysis as well.

```bash
magma --bfile g1000_eur --seed 42 --pval ../../iPSYCH-PGC_ASD_Nov2017 N=46351 \
    --gene-annot annot_WNH.genes.annot --out genes_ASD_WNH;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_ASD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_ASD_${r};
done;
```

### BD

I only found GWAS for multi-ethnic population only, which works for the main analysis.

```bash
# only grab rows we'll need
tail -n +73 ../../pgc-bip2021-all.vcf.tsv > ../../BD.txt;
# remove number sign manually, change SNP and P header

magma --bfile g1000_BW --seed 42 --pval ../../BD.txt N=51710 \
    --gene-annot annot_BW.genes.annot --out genes_BD_BW;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_BD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_BD_${r};
done;
```

### SCZ

It has EUR and multi-ethnic GWAS results. For the main analysis we only use the
latter.

```bash
magma --bfile g1000_BW --seed 42 \
    --pval ../../PGC3_SCZ_wave3_public.clumped.v2.tsv \
    N=161405 --gene-annot annot_BW.genes.annot --out genes_SCZ_BW;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_SCZ_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_SCZ_${r};
done;
```

### MDD

I only found GWAS for multi-ethnic population only, which works for the main analysis.

```bash
# changed SNP in header
magma --bfile g1000_BW --seed 42 \
    --pval ../../PGC_UKB_depression_genome-wide.txt N=500199 \
    --gene-annot annot_BW.genes.annot --out genes_MDD_BW;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_MDD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_MDD_${r};
done;
```

### OCD

I only found GWAS for multi-ethnic population only, which works for the main
analysis.

```bash
magma --bfile g1000_BW --seed 42 --pval ../../ocd_aug2017 N=9725 \
    --gene-annot annot_BW.genes.annot --out genes_OCD_BW;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_OCD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_OCD_${r};
done;
```

### PTSD

It has EUR and multi-ethnic GWAS results. For the main analysis we only use the
latter.

```bash
magma --bfile g1000_BW --seed 42 --pval ../../pts_all_freeze2_overall.results \
    N=206655 --gene-annot annot_BW.genes.annot --out genes_PTSD_BW;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_PTSD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_PTSD_${r};
done;
```

### Tourette

I only found GWAS for WNH population only. So I'll have to use the _eur BIM file
for the main analysis as well.

```bash
magma --bfile g1000_eur --seed 42 --pval ../../TS_Oct2018 N=14307 \
    --gene-annot annot_WNH.genes.annot --out genes_TS_WNH;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_TS_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_TS_${r};
done;
```

### Alcohol Use Disorder

I only found GWAS for WNH population only. So I'll have to use the _eur BIM file
for the main analysis as well.

```bash
# grab only rsids
head -n 1 ../../AUDIT_UKB_2018_AJP.txt > ../../AUD.txt;
grep rs ../../AUDIT_UKB_2018_AJP.txt >> ../../AUD.txt

# manually change SNP and p_T header (AUDIT total score)

magma --bfile g1000_eur --seed 42 --pval ../../AUD.txt N=121604 \
    --gene-annot annot_WNH.genes.annot --out genes_AUD_WNH;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_AUD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_AUD_${r};
done;
```

### Collecting results

Let's collect all results in R and save them in a matrix for plotting:

```r
g = 'main'
mydir = '~/data/post_mortem_adhd/results/MAGMA/'
dis = c('ADHD', 'ASD', 'OCD', 'BD', 'SCZ', 'TS', 'MDD', 'PTSD', 'AUD')
res = data.frame()
for (d in dis) {
    for (r in c('ACC', 'Caudate')) {
        fname = sprintf('%s/MAGMA_%s_gc_dge_%s_%s.gsa.out',
                        mydir, g, d, r)
        a = read.table(fname, header=1)
        b = a[1, 3:7]
        b$DISORDER = d
        b$POP = g
        b$REGION = r
        res = rbind(res, b)
    }
}
saveRDS(res, file=sprintf('~/data/post_mortem_adhd/results/MAGMA_%s_res.rds', g))
```


## Disorder correlation

```r
g = 'main'

do_boot_corrs = function(both_res, log2FC_col, method) {
    corrs = c()
    nperms = 10000
    set.seed(42)
    options(warn=-1)  # remove annoying spearman warnings
    for (p in 1:nperms) {
        idx = sample(nrow(both_res), replace = T)
        corrs = c(corrs, cor.test(abs(both_res[idx, 'log2FoldChange']),
                                  abs(both_res[idx, log2FC_col]),
                                  method=method)$estimate)
    }
    return(corrs)
}

out_fname = sprintf('~/data/post_mortem_adhd/results/disorders_corrs_%s.rds', g)
load(sprintf('~/data/post_mortem_adhd/results/%s_DGE.RData', g))

library(DESeq2)
meta = readRDS('~/data/post_mortem_adhd/aad6469_Gandal_SM_Data-Table-S1_micro.rds')

met = 'spearman'
dge = as.data.frame(results(dds.ACC, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by='ensembl_gene_id', all.x=F, all.y=F)

corrs = list()
disorders = c('ASD', 'SCZ', 'BD', 'MDD', 'AAD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('%s.beta_log2FC', d), met)
}
all_corrs = c()
for (d in disorders) {
    cat(d, '\n')
    junk = data.frame(corr=corrs[[d]])
    junk$region = 'ACC'
    junk$disorder = d
    junk$gene_overlap = nrow(both_res)
    junk$source = 'Gandal_micro'
    all_corrs = rbind(all_corrs, junk)
}

dge = as.data.frame(results(dds.Caudate, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by='ensembl_gene_id', all.x=F, all.y=F)
corrs = list()
disorders = c('ASD', 'SCZ', 'BD', 'MDD', 'AAD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('%s.beta_log2FC', d), met)
}
for (d in disorders) {
    junk = data.frame(corr=corrs[[d]])
    junk$region = 'Caudate'
    junk$disorder = d
    junk$gene_overlap = nrow(both_res)
    junk$source = 'Gandal_micro'
    all_corrs = rbind(all_corrs, junk)
}

library(gdata)
meta = read.xls('~/data/post_mortem_adhd/aad6469_Gandal_SM_Data-Table-S1.xlsx',
                'RNAseq SCZ&BD MetaAnalysis DGE')
dge = as.data.frame(results(dds.ACC, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='X', all.x=F, all.y=F)
corrs = list()
disorders = c('SCZ', 'BD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('%s.logFC', d), met)
}
for (d in disorders) {
    junk = data.frame(corr=corrs[[d]])
    junk$region = 'ACC'
    junk$disorder = d
    junk$gene_overlap = nrow(both_res)
    junk$source = 'Gandal_RNAseq'
    all_corrs = rbind(all_corrs, junk)
}

dge = as.data.frame(results(dds.Caudate, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='X', all.x=F, all.y=F)
corrs = list()
disorders = c('SCZ', 'BD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('%s.logFC', d), met)
}
for (d in disorders) {
    junk = data.frame(corr=corrs[[d]])
    junk$region = 'Caudate'
    junk$disorder = d
    junk$gene_overlap = nrow(both_res)
    junk$source = 'Gandal_RNAseq'
    all_corrs = rbind(all_corrs, junk)
}

meta = read.xls('~/data/post_mortem_adhd/aad6469_Gandal_SM_Data-Table-S1.xlsx',
                'RNAseq ASD-pancortical DGE')
dge = as.data.frame(results(dds.ACC, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='X', all.x=F, all.y=F)
corrs = list()
d = 'ASD'
junk = data.frame(corr=do_boot_corrs(both_res, 'Frontal.logFC', met))
junk$region = 'ACC'
junk$disorder = d
junk$gene_overlap = nrow(both_res)
junk$source = 'Gandal_RNAseq'
all_corrs = rbind(all_corrs, junk)

dge = as.data.frame(results(dds.Caudate, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='X', all.x=F, all.y=F)
corrs = list()
d = 'ASD'
junk = data.frame(corr=do_boot_corrs(both_res, 'Frontal.logFC', met))
junk$region = 'Caudate'
junk$disorder = d
junk$gene_overlap = nrow(both_res)
junk$source = 'Gandal_RNAseq'
all_corrs = rbind(all_corrs, junk)

# moving on to other papers: Akula
meta = readRDS('~/data/post_mortem_adhd/ACC_other_disorders.rds')
dge = as.data.frame(results(dds.ACC, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='Ensemble.gene.ID',
                 all.x=F, all.y=F)
corrs = list()
disorders = c('BD', 'SCZ', 'MDD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('log2FoldChange.%s', d), met)
}
for (d in disorders) {
    junk = data.frame(corr=corrs[[d]])
    junk$region = 'ACC'
    junk$disorder = d
    junk$gene_overlap = nrow(both_res)
    junk$source = 'Akula'
    all_corrs = rbind(all_corrs, junk)
}

dge = as.data.frame(results(dds.Caudate, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
mart = readRDS('~/data/post_mortem_adhd/mart_rnaseq.rds')
d = 'SCZ'
dge = merge(dge, mart, by='ensembl_gene_id', all.x=T, all.y=F)
meta = read.xls('~/data/post_mortem_adhd/caudate_others.xlsx', d)
meta$gencodeID = substr(meta$gencodeID, 1, 15)
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='gencodeID',
                 all.x=T, all.y=F)
colnames(both_res)[ncol(both_res)] = 'log2FC.SCZ'
junk = data.frame(corr=do_boot_corrs(both_res, sprintf('log2FC.%s', d), met))
junk$region = 'Caudate'
junk$disorder = d
junk$gene_overlap = nrow(both_res)
junk$source = 'Benjamin'
all_corrs = rbind(all_corrs, junk)

d = 'BD'
meta = read.xls('~/data/post_mortem_adhd/caudate_others.xlsx', d)
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='gencodeID',
                 all.x=T, all.y=F)
colnames(both_res)[ncol(both_res)] = 'log2FC.BD'
junk = data.frame(corr=do_boot_corrs(both_res, sprintf('log2FC.%s', d), met))
junk$region = 'Caudate'
junk$disorder = d
junk$gene_overlap = nrow(both_res)
junk$source = 'Pacifico'
all_corrs = rbind(all_corrs, junk)

d = 'OCD'
meta = read.xls('~/data/post_mortem_adhd/caudate_others.xlsx', d)
both_res = merge(dge, meta, by='hgnc_symbol', all.x=T, all.y=F)
colnames(both_res)[ncol(both_res)] = 'log2FC.OCD'
junk = data.frame(corr=do_boot_corrs(both_res, sprintf('log2FC.%s', d), met))
junk$region = 'Caudate'
junk$disorder = d
junk$gene_overlap = nrow(both_res)
junk$source = 'Piantadosi'
all_corrs = rbind(all_corrs, junk)

# last 2 ASD papers
dge = as.data.frame(results(dds.ACC, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
meta = read.xls('~/data/post_mortem_adhd/ASD_only.xlsx', 'Wright')
both_res = merge(dge, meta, by='ensembl_gene_id', all.x=T, all.y=F)
d = 'ASD'
junk = data.frame(corr=do_boot_corrs(both_res, 'log2FC', met))
junk$region = 'ACC'
junk$disorder = d
junk$gene_overlap = nrow(both_res)
junk$source = 'Wright_DLPFC'
all_corrs = rbind(all_corrs, junk)

meta = read.xls('~/data/post_mortem_adhd/ASD_only.xlsx', 'Neelroop')
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='ENSEMBL.ID',
                 all.x=T, all.y=F)
junk = data.frame(corr=do_boot_corrs(both_res, 'log2.FC..ASD.vs.CTL', met))
junk$region = 'ACC'
junk$disorder = d
junk$gene_overlap = nrow(both_res)
junk$source = 'Neelroop_FrontalTemporal'
all_corrs = rbind(all_corrs, junk)

saveRDS(all_corrs, file=out_fname)
```

To calculate p-values, we replace the bootstrap function and change the file
name, just to not have to repeat the code above:

```r
do_boot_corrs = function(both_res, log2FC_col, method) {
    corrs = c()
    nperms = 10000
    set.seed(42)
    options(warn=-1)  # remove annoying spearman warnings
    for (p in 1:nperms) {
        idx = sample(nrow(both_res), replace = F)
        corrs = c(corrs, cor.test(abs(both_res[, 'log2FoldChange']),
                                  abs(both_res[idx, log2FC_col]),
                                  method=method)$estimate)
    }
    return(corrs)
}
out_fname = sprintf('~/data/post_mortem_adhd/results/disorders_corrs_null_%s.rds',
                    g)
```


# Robustness analysis

## White non-Hispanic only analysis

### DGE

Same code as main analysis, except for: 

```r
g = 'WNH'
dds.ACC = pca_DGE_clean('ACC', fm_str, varvec, wnhOnly=T)
dds.Caudate = pca_DGE_clean('Caudate', fm_str, varvec, wnhOnly=T)
```

and keep the change to g to export.


### GSEA

We don't need to re-create any of the GMT files. All we need is to change g in
the main code.

### MAGMA

For MAGMA we need to use the appropriate population.

#### Generic files

```bash
cd ~/data/post_mortem_adhd/results/MAGMA;
module load MAGMA/1.09a;

# annotation step
magma --annotate --seed 42 --snp-loc g1000_eur.bim \
    --gene-loc /usr/local/apps/MAGMA/gene_location/NCBI37.3/NCBI37.3.gene.loc \
    --out annot_WNH;
```

To create the gene covariance files, the only changes from above are:

```r
g = 'WNH'
```

#### ADHD

It has EUR and all GWAS results. Here we only need EUR:

```bash
g=WNH;

magma --bfile g1000_eur --seed 42 --pval ../../adhd_eur_jun2017 N=53293 \
    --gene-annot annot_WNH.genes.annot --out genes_ADHD_WNH;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_ADHD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_ADHD_${r};
done;
```

#### ASD

I only found GWAS for WNH population only, which I already used in the
construction step for the main analysis:

```bash
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_ASD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_ASD_${r};
done;
```

#### BD

I only found GWAS for multi-ethnic population only, which I already used in the
construction step for the main analysis:

```bash
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_BD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_BD_${r};
done;
```

#### SCZ

Need to prepare the file for WNH population:

```bash
# only use the rs SNPs
head -n 1 ../../daner_natgen_pgc_eur > ../../SCZ.txt;
grep rs ../../daner_natgen_pgc_eur >> ../../SCZ.txt;
magma --bfile g1000_eur --seed 42 --pval ../../SCZ.txt N=147521 \
    --gene-annot annot_WNH.genes.annot --out genes_SCZ_WNH
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_SCZ_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_SCZ_${r};
done;
```

#### MDD

I only found GWAS for multi-ethnic population only, which I already used in the
construction step for the main analysis:

```bash
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_MDD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_MDD_${r};
done;
```

####  OCD

I only found GWAS for multi-ethnic population only, which I already used in the
construction step for the main analysis:

```bash
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_OCD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_OCD_${r};
done;
```

#### PTSD

Still need to prepare the out file for EUR:

```bash
magma --bfile g1000_eur --seed 42 --pval ../../pts_eur_freeze2_overall.results \
    N=174659 --gene-annot annot_WNH.genes.annot --out genes_PTSD_WNH;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_PTSD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_PTSD_${r};
done;
```

#### Tourette

I only found GWAS for WNH population only, which I already used in the
construction step for the main analysis:

```bash
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_TS_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_TS_${r};
done;
```

#### Alcohol Use Disorder

I only found GWAS for WNH population only, which I already used in the
construction step for the main analysis:

```bash
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_AUD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_AUD_${r};
done;
```

#### Collecting results

Same code as above, but change g.

### Disorder correlation

Same code as above, except for changing g.

## Not using comorbidities and substance abuse as covariates

### DGE

Same code as main analysis, except for: 

```r
g = 'lessCov'
fm_str = '~ RINe + C1 + BBB2 + Diagnosis'
dds.ACC = pca_DGE_clean('ACC', fm_str, varvec, wnhOnly=F)
dds.Caudate = pca_DGE_clean('Caudate', fm_str, varvec, wnhOnly=F)
```

and keep the change to g to export.

### GSEA

We don't need to re-create any of the GMT files. All we need is to change g in
the main code.

### MAGMA

For MAGMA we will use the same population as the one for the main analysis:

#### Generic files

```r
g = 'lessCov'
```

#### ADHD

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_ADHD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_ADHD_${r};
done;
```

#### ASD

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_ASD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_ASD_${r};
done;
```

#### BD

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_BD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_BD_${r};
done;
```

#### SCZ

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_SCZ_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_SCZ_${r};
done;
```

#### MDD

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_MDD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_MDD_${r};
done;
```

####  OCD

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_OCD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_OCD_${r};
done;
```

#### PTSD

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_PTSD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_PTSD_${r};
done;
```

#### Tourette

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_TS_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_TS_${r};
done;
```

#### Alcohol Use Disorder

```bash
g=lessCov;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_AUD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_AUD_${r};
done;
```

#### Collecting results

Same main code above, but change g to lessCov.

### Disorder correlation

Same main code above, but change g to lessCov.


## Removing samples of subjects with MDD

### DGE

Same code as main analysis, except for: 

```r
g = 'noMDD'
fm_str = '~ RINe + C1 + BBB2 + comorbid_group + SUB2 + Diagnosis'
# list values to keep for that particular variable
varvec = list(comorbid_update=c('nil', 'adjustment disorder ',
                                'BPAD_NOS; Impulse Control Disorder', 'ASD',
                                "dysthmia"))
dds.ACC = pca_DGE_clean('ACC', fm_str, varvec, wnhOnly=F)
dds.Caudate = pca_DGE_clean('Caudate', fm_str, varvec, wnhOnly=F)
```

and keep the change to g to export.

### GSEA

We don't need to re-create any of the GMT files. All we need is to change g in
the main code.

### MAGMA

For MAGMA we will use the same BW population as the one for the main analysis:

#### Generic files

```r
g = 'noMDD'
```

#### ADHD

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_ADHD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_ADHD_${r};
done;
```

#### ASD

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_ASD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_ASD_${r};
done;
```

#### BD

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_BD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_BD_${r};
done;
```

#### SCZ

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_SCZ_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_SCZ_${r};
done;
```

#### MDD

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_MDD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_MDD_${r};
done;
```

####  OCD

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_OCD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_OCD_${r};
done;
```

#### PTSD

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_PTSD_BW.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_PTSD_${r};
done;
```

#### Tourette

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_TS_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_TS_${r};
done;
```

#### Alcohol Use Disorder

```bash
g=noMDD;
for r in 'ACC' 'Caudate'; do
    magma --seed 42 --gene-results genes_AUD_WNH.genes.raw \
        --gene-covar ../MAGMA_${g}_dge_${r}.tab \
        --out MAGMA_${g}_gc_dge_AUD_${r};
done;
```

#### Collecting results

Same main code above, but change g to noMDD.

### Disorder correlation

Same main code above, but change g to noMDD.


## Leave-one-out

Let's save only the final table for each subject left out:

```r
library(DESeq2)
mydir = '~/data/post_mortem_adhd/'
myregion = 'Caudate' # ACC or Caudate
mart = readRDS(sprintf('%s/mart_rnaseq.rds', mydir))

library(GenomicFeatures)
txdb <- loadDb(sprintf('%s/Homo_sapies.GRCh38.97.sqlite', mydir))
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
bt = read.csv(sprintf('%s/Homo_sapiens.GRCh38.97_biotypes.csv', mydir))
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]

set.seed(42)
fm_str = '~ RINe + C1 + BBB2 + comorbid_group + SUB2 + Diagnosis'
nsubjs = ifelse(myregion == 'ACC', 53, 56)
for (s in 1:nsubjs) {
    fname = sprintf('%s/results/LOO/DGE_%s_main_annot_loo%02d.csv', mydir,
                    myregion, s)
    dds = pca_DGE_clean(myregion, fm_str, NA, loo=s)
    res = results(dds, name = "Diagnosis_Case_vs_Control", alpha=.05)
    
    df = as.data.frame(res)
    colnames(df)[ncol(df)] = 'padj.FDR'
    # grabbing HUGO IDs
    df$GENEID = substr(rownames(df), 1, 15)
    df2 = merge(df, mart, sort=F,
                by.x='GENEID', by.y='ensembl_gene_id', all.x=T, all.y=F)
    df2 = merge(df2, bt_slim, sort=F,
                by.x='GENEID', by.y='gene_id', all.x=T, all.y=F)
    df2 = df2[order(df2$pvalue), ]
    
    write.csv(df2, row.names=F, file=fname)
}
```

Now we check the rank and other stats for each of the significant DEGs:

```r
library(DESeq2)
mydir = '~/data/post_mortem_adhd/results/'
load(sprintf('%s/main_DGE.RData', mydir))
for (r in c('ACC', 'Caudate')) {
    fname = sprintf('%s/DGE_%s_main_annot.csv', mydir, r)
    res = read.csv(fname)
    top_genes = res[round(res$padj.FDR, 2) <= .05, 'GENEID']
    hgnc = res[round(res$padj.FDR, 2) <= .05, 'hgnc_symbol']
    all_res = data.frame(GENEID=top_genes, HUGOID=hgnc, ES=0,
                         rank=1:length(top_genes), p_boot=0, q_boot=0)

    # following this: https://www.biostars.org/p/140976/ for standardized ES
    res = res[round(res$padj.FDR, 2) <= .05, ]
    all_res$ES = res$log2FoldChange
    res_str = sprintf('dds = dds.%s', r)
    eval(parse(text=res_str))
    tmp = data.frame(dispersion=mcols(dds)$dispGeneEst)
    tmp$GENEID = substring(rownames(mcols(dds)), 1, 15)
    res2 = merge(res, tmp, by='GENEID', sort=F)
    res2$stdES = res2$log2FoldChange / sqrt(1/res2$baseMean + res2$dispersion)
    all_res$stdES = res2$stdES

    # using this for CIs: https://support.bioconductor.org/p/80725/
    res2$lCI = res2$log2FoldChange + qnorm(.025)*res2$lfcSE
    res2$uCI = res2$log2FoldChange + qnorm(.975)*res2$lfcSE
    all_res$lCI = res2$lCI
    all_res$uCI = res2$uCI

    all_res = all_res[order(all_res$GENEID), ]

    # collecting leave-one-out results
    boot_files = list.files(path=sprintf('%s/LOO/', mydir),
                            pattern=sprintf('^DGE_%s_main_annot_loo', r))
    boot_ranks = matrix(nrow=length(top_genes), ncol=length(boot_files))
    for (f in 1:length(boot_files)) {
        fname = sprintf('%s/LOO/%s', mydir, boot_files[f])
        res = read.csv(fname)
        for (g in 1:length(top_genes)) {
            idx = res$GENEID == top_genes[g]
            if (any(idx)) {
                boot_ranks[g, f] = which(idx)
            }
        }        
        res = res[res$GENEID %in% top_genes, ]
        res = res[order(res$GENEID),]
        if (all(all_res$GENEID == res$GENEID)) {
            idx = which(round(res$pvalue, 2) <= .05)
            all_res[idx, 'p_boot'] = all_res[idx, 'p_boot'] + 1
            idx = which(round(res$padj.FDR, 2) <= .05)
            all_res[idx, 'q_boot'] = all_res[idx, 'q_boot'] + 1
        } else {
            cat('LOO', f, 'does not have all genes\n')
        }
    }
    all_res$med_rank = apply(boot_ranks, 1, median, na.rm=T)

    for (i in 1:nrow(all_res)) {
        all_res[i, 'Iterations p <= .05'] = sprintf('%.0f / %.0f',
                                                   all_res[i, 'p_boot'],
                                                   length(boot_files))
        all_res[i, 'Iterations q <= .05'] = sprintf('%.0f / %.0f',
                                                   all_res[i, 'q_boot'],
                                                   length(boot_files))
        all_res[i, 'Effect size (95%% CI)'] = sprintf('%.2f (%.2f; %.2f)',
                                                     all_res[i, 'ES'],
                                                     all_res[i, 'lCI'],
                                                     all_res[i, 'uCI'])
    }

    colnames(all_res)[c(4, 10)] = c('Rank', 'Median LOO Rank')
    all_res = all_res[order(all_res$Rank), ]
    out_fname = sprintf('%s/DGE_%s_main_LOO.csv', mydir, r)
    write.csv(all_res[, c(1, 2, 4, 10:13)], file=out_fname, row.names=F, quote=F)
    print(all_res)
}
```

## Session info

Just recording the version of all packages we used in the analysis:

```
r$> sessionInfo()                                                                                                                  
R version 4.1.3 (2022-03-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/local/intel/compilers_and_libraries_2020.2.254/linux/mkl/lib/intel64_lin/libmkl_rt.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dplyr_1.0.8                 org.Hs.eg.db_3.14.0         WebGestaltR_0.4.4           ABAEnrichment_1.24.0       
 [5] ActivePathways_1.1.0        edgeR_3.36.0                limma_3.50.1                DESeq2_1.34.0              
 [9] SummarizedExperiment_1.24.0 MatrixGenerics_1.6.0        matrixStats_0.61.0          GenomicFeatures_1.46.5     
[13] AnnotationDbi_1.56.2        Biobase_2.54.0              GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
[17] IRanges_2.28.0              S4Vectors_0.32.4            BiocGenerics_0.40.0         gdata_2.18.0               

loaded via a namespace (and not attached):
 [1] bitops_1.0-7             bit64_4.0.5              doParallel_1.0.17        filelock_1.0.2           RColorBrewer_1.1-3      
 [6] progress_1.2.2           httr_1.4.2               doRNG_1.8.2              tools_4.1.3              utf8_1.2.2              
[11] R6_2.5.1                 KernSmooth_2.23-20       DBI_1.1.2                colorspace_2.0-3         tidyselect_1.1.2        
[16] prettyunits_1.1.1        bit_4.0.4                curl_4.3.2               compiler_4.1.3           ABAData_1.24.0          
[21] cli_3.2.0                xml2_1.3.3               DelayedArray_0.20.0      rtracklayer_1.54.0       caTools_1.18.2          
[26] scales_1.1.1             readr_2.1.2              genefilter_1.76.0        rappdirs_0.3.3           systemfonts_1.0.4       
[31] apcluster_1.4.9          stringr_1.4.0            digest_0.6.29            Rsamtools_2.10.0         svglite_2.1.0           
[36] XVector_0.34.0           pkgconfig_2.0.3          dbplyr_2.1.1             fastmap_1.1.0            rlang_1.0.2             
[41] RSQLite_2.2.12           BiocIO_1.4.0             generics_0.1.2           jsonlite_1.8.0           BiocParallel_1.28.3     
[46] gtools_3.9.2             RCurl_1.98-1.6           magrittr_2.0.3           GenomeInfoDbData_1.2.7   Matrix_1.4-1            
[51] Rcpp_1.0.8.3             munsell_0.5.0            fansi_1.0.3              lifecycle_1.0.1          whisker_0.4             
[56] stringi_1.7.6            yaml_2.3.5               zlibbioc_1.40.0          gplots_3.1.1             BiocFileCache_2.2.1     
[61] grid_4.1.3               blob_1.2.2               parallel_4.1.3           crayon_1.5.1             lattice_0.20-45         
[66] Biostrings_2.62.0        splines_4.1.3            annotate_1.72.0          hms_1.1.1                KEGGREST_1.34.0         
[71] locfit_1.5-9.5           pillar_1.7.0             igraph_1.3.0             rjson_0.2.21             rngtools_1.5.2          
[76] codetools_0.2-18         geneplotter_1.72.0       biomaRt_2.50.3           XML_3.99-0.9             glue_1.6.2              
[81] data.table_1.14.2        tzdb_0.3.0               foreach_1.5.2            png_0.1-7                vctrs_0.4.0             
[86] gtable_0.3.0             purrr_0.3.4              assertthat_0.2.1         cachem_1.0.6             ggplot2_3.3.5           
[91] xtable_1.8-4             restfulr_0.0.13          survival_3.3-1           tibble_3.1.6             iterators_1.0.14        
[96] GenomicAlignments_1.30.0 memoise_2.0.1            ellipsis_0.3.2          

```

# Exporting results

Let's export all our results into nice Excel tables, with a glossary:

```r
g = 'noMDD'

res_dir = '~/data/post_mortem_adhd/results/'
out_fname = sprintf('%s/eFile_results_%s.xlsx', res_dir, g)
library(openxlsx)
wb = createWorkbook()
sHeader <- createStyle(textDecoration = c("BOLD"))
for (r in c('ACC', 'Caudate')) {
    in_fname = sprintf('%s/DGE_%s_%s_annot.csv', res_dir, r, g)
    sname = sprintf('DGE_%s', r)
    addWorksheet(wb, sname)
    res = read.csv(in_fname)
    writeData(wb, sheet=sname, res, rowNames=F)
    addStyle(wb, sheet = sname, sHeader, rows=1, cols=1:ncol(res))
    
    DBs = c('Biological_Process', 'Cellular_Component', 'Molecular_Function')
    for (db in DBs) {
        in_fname = sprintf('%s/GSEA_%s_%s_geneontology_%s_noRedundant_10K.csv',
                           res_dir, g, r, db)
        sname = sprintf('GSEA_%s_%s', r, db)
        addWorksheet(wb, sname)
        res = read.csv(in_fname)
        res$plotPath = NULL
        writeData(wb, sheet=sname, res, rowNames=F)
        addStyle(wb, sheet = sname, sHeader, rows=1, cols=1:ncol(res))
    }
    db = 'developmental'
    in_fname = sprintf('%s/GSEA_%s_%s_%s_10K.csv',
                       res_dir, g, r, db)
    sname = sprintf('GSEA_%s_%s', r, db)
    addWorksheet(wb, sname)
    res = read.csv(in_fname)
    res$plotPath = NULL
    writeData(wb, sheet=sname, res, rowNames=F)
    addStyle(wb, sheet = sname, sHeader, rows=1, cols=1:ncol(res))
}
saveWorkbook(wb, file=out_fname)
```

And we also add a glossary for the different columns:

```r
res = c()

res = rbind(res, c('DGE', 'GENEID', 'GENE Ensemble ID'))
res = rbind(res, c('DGE', 'baseMean', 'mean of normalized counts for all samples'))
res = rbind(res, c('DGE', 'log2FoldChange', 'log2 fold change (MLE): ADHD vs Unaffected'))
res = rbind(res, c('DGE', 'lfcSE', 'standard error: ADHD vs Unaffected'))
res = rbind(res, c('DGE', 'stat', 'Wald statistic: ADHD vs Unaffected'))
res = rbind(res, c('DGE', 'pvalue', 'Wald test p-value: ADHD vs Unaffected'))
res = rbind(res, c('DGE', 'padj.FDR', 'BH adjusted p-values'))
res = rbind(res, c('DGE', 'hgnc_symbol', 'HUGO gene symbol'))
res = rbind(res, c('DGE', 'chromosome_name', 'Chromosome where gene is located'))
res = rbind(res, c('DGE', 'gene_biotype', 'gene type'))

res = rbind(res, c('GSEA', 'geneSet', 'ID of the gene set'))
res = rbind(res, c('GSEA', 'description', 'Desription of the gene set'))
res = rbind(res, c('GSEA', 'link', 'Link to the data source'))
res = rbind(res, c('GSEA', 'enrichmentScore', 'The maximum running sum of scores for the ranked list'))
res = rbind(res, c('GSEA', 'normalizedEnrichmentScore', 'Enrichment score normalized against the average enrichment score of all permutations'))
res = rbind(res, c('GSEA', 'pValue', 'Nominal p-value. The statistical significance of the enrichment score. The nominal p value is not adjusted for gene set size or multiple hypothesis testing; therefore, it is of limited use in comparing gene sets (see https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm)'))
res = rbind(res, c('GSEA', 'FDR', 'False discovery rate; that is, the estimated probability that the normalized enrichment score represents a false positive finding (see https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm).'))
res = rbind(res, c('GSEA', 'size', 'The number of genes in the set after filtering by minNum and maxNum.'))
res = rbind(res, c('GSEA', 'leadingEdgeNum', 'Number of genes in the leading edge subset. Those are genes that appear in the ranked list at or before the point at which the running sum reaches its maximum deviation from zero. The leading-edge subset can be interpreted as the core that accounts for the gene set’s enrichment signal.'))
res = rbind(res, c('GSEA', 'leadingEdgeId', 'Genes in the leading edge subset in Entrez geneID'))
res = rbind(res, c('GSEA', 'userId', 'Genes in the leading edge subset in Ensemble ID'))

colnames(res) = c('tab', 'column', 'description')

wb = createWorkbook()
sHeader <- createStyle(textDecoration = c("BOLD"))
sname = 'main'
addWorksheet(wb, sname)
writeData(wb, sheet=sname, res, rowNames=F)
addStyle(wb, sheet = sname, sHeader, rows=1, cols=1:ncol(res))
out_fname = sprintf('%s/eFile_glossary.xlsx', res_dir)
saveWorkbook(wb, file=out_fname)
```

# Figures

## Figure 1

### Volcano plots

```r
lsize = 12
tsize = 9
msize = 16

library(ggpubr)
library(EnhancedVolcano)
FCcutoff = 1.0
pCutoff = .05

myplots = list()
res = read.csv('~/data/post_mortem_adhd/results/DGE_ACC_main_annot.csv')
res = res[order(res$pvalue), ]
sigPs = sum(res$padj.FDR <= pCutoff, na.rm=T)
ps = -log10(res$pvalue)
nomPcutoff = ps[sigPs + 1] + (ps[sigPs] - ps[sigPs + 1]) / 2
nomPcutoff = 10 ** (-nomPcutoff)
ymax = ceiling(max(-log10(res$pvalue), na.rm=T))
xmin = floor(min(res$log2FoldChange, na.rm=T))
xmax = ceiling(max(res$log2FoldChange, na.rm=T))
p = EnhancedVolcano(data.frame(res),
                    x = 'log2FoldChange', lab=NA,
                    y = 'pvalue', xlab = bquote(~log[2]~ 'fold change'),
                    title = 'ACC', titleLabSize=msize, 
                    ylab = bquote(~-log[10]~italic((p))),
                    ylim = c(0, ymax),
                    xlim = c(xmin, xmax),
                    pCutoff = nomPcutoff, FCcutoff = FCcutoff, pointSize = 1.0,
                    subtitle=NULL,
                    axisLabSize = lsize,
                    caption = NULL, legendPosition = 'none',
                    col = c("grey30", "grey30", "#DC3220", "#DC3220"))
myplots[[1]] = p

res = read.csv('~/data/post_mortem_adhd/results/DGE_Caudate_main_annot.csv')
sigPs = sum(res$padj.FDR <= pCutoff, na.rm=T)
ps = -log10(res$pvalue)
nomPcutoff = ps[sigPs + 1] + (ps[sigPs] - ps[sigPs + 1]) / 2
nomPcutoff = 10 ** (-nomPcutoff)
p = EnhancedVolcano(data.frame(res),
                    x = 'log2FoldChange', lab=NA,
                    y = 'pvalue', xlab = bquote(~log[2]~ 'fold change'),
                    title = 'Caudate', titleLabSize=msize, 
                    ylim = c(0, ymax),
                    xlim = c(xmin, xmax),
                    pCutoff = nomPcutoff, FCcutoff = FCcutoff, pointSize = 1.0,
                    subtitle=NULL,
                    axisLabSize = lsize,
                    caption = NULL, legendPosition = 'none',
                    col = c("grey30", "grey30", "#DC3220", "#DC3220"))
myplots[[2]] = p + theme(axis.title.y = element_blank())

volcanoes = ggarrange(plotlist=myplots)
```

### Manhattan plots

```r
# using a heavily eddited version of the script found in https://github-wiki-see.page/m/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function
# because it allows to use ggrepel and avoid overlapping labels

library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

gg.manhattan <- function(df, threshold, hlight, ylims, cols, title, xlab, ylab){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(TXSTART)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, TXSTART) %>%
    mutate( BPcum=TXSTART+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(hgnc_symbol %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(pvalue < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  p = ggplot(df.tmp, aes(x=BPcum, y=-log10(pvalue))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1) +
    scale_color_manual(values = rep(cols, 22 )) +

    # custom X axis:
    scale_x_continuous( label = c('1', '2', '3', '4', '5', '6', '7', '8', '9',
                                  '10', '11', '12', '13', '14', '15', '16', '',
                                  '18', '', '20', '', '22'),
                                  breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = xlab, y=ylab) +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(threshold), color='#DC3220', linetype='dashed') +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(hgnc_symbol), alpha=0.2), size=3, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5, size = msize, face='bold'),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_text(size = lsize),
      axis.text.y = element_text(size = tsize),
      axis.title.x = element_text(size = lsize),
      axis.text.x = element_text(size = tsize)
    )
    return(p)
}

library(AnnotationDbi)
library(DESeq2)

load('~/data/post_mortem_adhd/results/main_DGE.RData')
dge = as.data.frame(results(dds.ACC, name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
mart = readRDS('~/data/post_mortem_adhd/mart_rnaseq.rds')
dge = merge(dge, mart, by='ensembl_gene_id', all.x=T, all.y=F)
txdb <- loadDb('~/data/post_mortem_adhd/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"),
               columns=c('GENEID','TXCHROM', 'TXSTART', 'TXEND'),
               "GENEID")
txdf = txdf[!duplicated(txdf$GENEID),] 
dge = merge(dge, txdf, by.x='ensembl_gene_id', by.y='GENEID', all.x=T,
            all.y=F)
dge$CHR = as.numeric(dge$chromosome_name)
dge = dge[!is.na(dge$CHR), ]
idx = dge$hgnc_symbol==''
dge[idx, 'hgnc_symbol'] = dge[idx, 'ensembl_gene_id']
# forcing it based on the FDR-adjusted p-values
p_thresh = 3e-05
p1 = gg.manhattan(dge, p_thresh, NA, c(0,10), c('black', 'grey'), "ACC", ' ',
                  bquote(~-log[10]~italic((p))))
p1 = p1 + theme(axis.title.x = element_blank())

dge = as.data.frame(results(dds.Caudate,
                    name = "Diagnosis_Case_vs_Control"))
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
dge = merge(dge, mart, by='ensembl_gene_id', all.x=T, all.y=F)
dge = merge(dge, txdf, by.x='ensembl_gene_id', by.y='GENEID', all.x=T,
            all.y=F)
dge$CHR = as.numeric(dge$chromosome_name)
dge = dge[!is.na(dge$CHR), ]
idx = dge$hgnc_symbol==''
dge[idx, 'hgnc_symbol'] = dge[idx, 'ensembl_gene_id']
p_thresh = 9.33e-07
p2 = gg.manhattan(dge, p_thresh, NA, c(0,10), c('black', 'grey'), "Caudate",
                  'Chromosome', bquote(~-log[10]~italic((p))))

manh = ggarrange(p1, p2, nrow=2, ncol=1)

fig = ggarrange(volcanoes, manh, nrow=2, ncol=1, labels='AUTO',
                heights=c(.8, 2.2))
# letter is 8.5 by 11 in
ggsave("~/data/post_mortem_adhd/figures/figure1.pdf", fig, width=7, height=9,
       units="in")
```

## Figure 2

### Molecular function plots (FDR q <= .05)

```r
tsize = 12
ysize = 10
xsize = 9
msize = 2
ntop = 10

df = read.csv('~/data/post_mortem_adhd/results/GSEA_main_Caudate_geneontology_Molecular_Function_noRedundant_10K.csv')
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore', 'FDR')]
df = df[round(df$FDR, 2) <= 0.05, ]
df[df$FDR == 0, 'FDR'] = 1e-4

df$description = factor(df$description,
                        levels=df$description[sort(df$FDR,
                                                   index.return=T,
                                                   decreasing=T)$ix])
color_me = c('fatty acid derivative binding', 'dynein light chain binding',
             'tau-protein kinase activity')
label_colors <- ifelse(levels(df$description) %in% color_me, "grey20", "#DC3220")
df$Behavior = ifelse(df$description %in% color_me, "notneuro", "neuro")

p <- ggplot(df, aes(y=-log10(FDR), x=description, fill=Behavior)) +
  geom_dotplot(dotsize=msize*.9, binaxis='y', stackdir='center') + coord_flip() +
  geom_hline(yintercept=-log10(.1), linetype="dashed",
                                color = "black", size=1) + theme(legend.position="bottom") +
    geom_hline(yintercept=-log10(.05), linetype="dotted",
                                color = "black", size=1) + theme(legend.position="bottom") +
    theme(axis.text.y = element_text(colour = label_colors, size=ysize),
          axis.title.y = element_blank(),
          plot.title = element_text(size = tsize),
          axis.text.x = element_text(size=xsize),
          axis.title.x = element_text(size=ysize)) +
    scale_fill_manual(values=c("#DC3220", 'grey60'))

library(ggpubr)
p1 = p + ggtitle('Caudate') + ylab(bquote(~-log[10]~italic((p[adjusted])))) +
    ylim(0, 4.25)

df = read.csv('~/data/post_mortem_adhd/results/GSEA_main_ACC_geneontology_Molecular_Function_noRedundant_10K.csv')
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore', 'FDR')]
df = df[round(df$FDR, 2) <= 0.05, ]
df[df$FDR == 0, 'FDR'] = 1e-4

df$description = factor(df$description,
                        levels=df$description[sort(df$FDR,
                                                   index.return=T,
                                                   decreasing=T)$ix])

color_me = c('neurotransmitter receptor activity',
             'serotonin receptor activity', 'GABA receptor activity')
label_colors <- ifelse(levels(df$description) %in% color_me, "#DC3220", "grey20")
df$Behavior = ifelse(df$description %in% color_me, "neuro", "notneuro")

p <- ggplot(df, aes(y=-log10(FDR), x=description, fill=Behavior, size=msize)) +
  geom_dotplot(dotsize=1.5*msize, binaxis='y', stackdir='center') + coord_flip() +
  geom_hline(yintercept=-log10(.1), linetype="dashed",
                                color = "black", size=1) + theme(legend.position="bottom") +
    geom_hline(yintercept=-log10(.05), linetype="dotted",
                                color = "black", size=1) + theme(legend.position="bottom") +
    theme(axis.text.y = element_text(colour = label_colors, size=ysize),
          axis.title.y = element_blank(),
          plot.title = element_text(size = tsize),
          axis.text.x = element_text(size=xsize),
          axis.title.x = element_text(size=ysize)) +
    scale_fill_manual(values=c("#DC3220", 'grey60'))

p2 = p + ggtitle('ACC') + ylab(bquote(~-log[10]~italic((p[adjusted])))) +
    ylim(0, 4.25)

mol = ggarrange(p1, p2, common.legend=F, ncol=2, nrow=1, legend='none')
```

### Top 10 Cellular components

```r
df = read.csv('~/data/post_mortem_adhd/results/GSEA_main_Caudate_geneontology_Cellular_Component_noRedundant_10K.csv')
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore', 'FDR')]
df = df[1:ntop, ]
df[df$FDR == 0, 'FDR'] = 1e-4

df$description = factor(df$description,
                        levels=df$description[sort(df$FDR,
                                                   index.return=T,
                                                   decreasing=T)$ix])
color_me = c('transporter complex', 'mitochondrial protein complex')
label_colors <- ifelse(levels(df$description) %in% color_me, "grey20", "#DC3220")
df$Behavior = ifelse(df$description %in% color_me, "notneuro", "neuro")

p <- ggplot(df, aes(y=-log10(FDR), x=description, fill=Behavior)) +
  geom_dotplot(dotsize=msize, binaxis='y', stackdir='center') + coord_flip() +
  geom_hline(yintercept=-log10(.1), linetype="dashed",
                                color = "black", size=1) + theme(legend.position="bottom") +
    geom_hline(yintercept=-log10(.05), linetype="dotted",
                                color = "black", size=1) + theme(legend.position="bottom") +
    theme(axis.text.y = element_text(colour = label_colors, size=ysize),
          axis.title.y = element_blank(),
          plot.title = element_text(size = tsize),
          axis.text.x = element_text(size=xsize),
          axis.title.x = element_text(size=ysize)) +
    scale_fill_discrete(drop = FALSE) +
    scale_fill_manual(values=c("#DC3220", 'grey60'))

library(ggpubr)
p1 = p + ggtitle('Caudate') + ylab(bquote(~-log[10]~italic((p[adjusted])))) +
     ylim(0, 5)

df = read.csv('~/data/post_mortem_adhd/results/GSEA_main_ACC_geneontology_Cellular_Component_noRedundant_10K.csv')
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore', 'FDR')]
df = df[1:ntop, ]
df[df$FDR == 0, 'FDR'] = 1e-4

df$description = factor(df$description,
                        levels=df$description[sort(df$FDR,
                                                   index.return=T,
                                                   decreasing=T)$ix])

color_me = c('nothing')
label_colors <- ifelse(levels(df$description) %in% color_me, "#DC3220", "grey20")
df$Behavior = ifelse(df$description %in% color_me, "neuro", "notneuro")
df$Behavior = factor(df$Behavior, levels=c('notneuro', 'neuro'))

p <- ggplot(df, aes(y=-log10(FDR), x=description, fill=Behavior, size=msize)) +
  geom_dotplot(dotsize=msize, binaxis='y', stackdir='center') + coord_flip() +
  geom_hline(yintercept=-log10(.1), linetype="dashed",
                                color = "black", size=1) + theme(legend.position="bottom") +
    geom_hline(yintercept=-log10(.05), linetype="dotted",
                                color = "black", size=1) + theme(legend.position="bottom") +
    theme(axis.text.y = element_text(colour = label_colors, size=ysize),
          axis.title.y = element_blank(),
          plot.title = element_text(size = tsize),
          axis.text.x = element_text(size=xsize),
          axis.title.x = element_text(size=ysize)) +
    scale_fill_discrete(drop = FALSE) +
    scale_fill_manual(values=c('grey60', "#DC3220"))

p2 = p + ggtitle('ACC') + ylab(bquote(~-log[10]~italic((p[adjusted])))) + 
     ylim(0, 5) 

cell = ggarrange(p1, p2, common.legend=T, ncol=2, nrow=1, legend='none')

fig = ggarrange(mol, cell, nrow=2, ncol=1, labels='AUTO',
                heights=c(.9, 1.1))
# letter is 8.5 by 11 in
ggsave("~/data/post_mortem_adhd/figures/figure2.pdf", fig, width=7.5, height=7,
       units="in")
```

## Figure 3

```r
library(corrplot)

r = 'ACC'
file_name = sprintf('~/data/post_mortem_adhd/results/GSEA_main_%s_developmental_10K.csv', r)
res = read.csv(file_name)
res = res[order(res$geneSet), c('link', 'normalizedEnrichmentScore', 'pValue')]
dev = res
dev$Region = r

df = matrix(nrow = 2, ncol = 6, dimnames=list(c('ACC', 'Caudate'),
                                              c('pandev', res$link[1:5])))
pvals = df
i = 1
for (j in 1:ncol(df)) {
    idx = dev$Region == rownames(df)[i] & dev$link == colnames(df)[j]
    if (dev[idx, 'pValue'] == 0) {
        dev[idx, 'pValue'] = 1e-4
    }
    df[i, j] = (sign(dev[idx, 'normalizedEnrichmentScore']) *
                -log10(dev[idx, 'pValue']))
    pvals[i, j] = dev[idx, 'pValue']
}

r = 'Caudate'
file_name = sprintf('~/data/post_mortem_adhd/results/GSEA_main_%s_developmental_10K.csv', r)
res = read.csv(file_name)
res = res[order(res$geneSet), c('link', 'normalizedEnrichmentScore', 'pValue')]
res$Region = r
dev = res

i = 2
for (j in 1:ncol(df)) {
    idx = dev$Region == rownames(df)[i] & dev$link == colnames(df)[j]
    if (dev[idx, 'pValue'] == 0) {
        dev[idx, 'pValue'] = 1e-4
    }
    df[i, j] = (sign(dev[idx, 'normalizedEnrichmentScore']) *
                -log10(dev[idx, 'pValue']))
    pvals[i, j] = dev[idx, 'pValue']
}
mylim = max(abs(df))
colnames(df)[1] = 'pan-developmental'

# i checked and it comes out in the correct order
pvals2 = matrix(p.adjust(pvals, method='fdr'), ncol=6, nrow=2)
colnames(pvals2) = colnames(df)
rownames(pvals2) = rownames(df)

# figuring out rows and columns to highlight
markme = which(round(pvals2, 2) <= .05, arr.ind = T) 

pdf(file="~/data/post_mortem_adhd/figures/figure3.pdf", width=5, height=3)
corrplot(df, is.corr=F, col.lim=c(-mylim, mylim), tl.col='black',
         cl.length=3, cl.align.text='l', cl.offset=.2,
         col=colorRampPalette(c("#DC3220","white","#005AB5"))(200),) -> p
# make thicker borders for each significant square
for (row in 1:nrow(markme)) {
    r = c(colnames(df)[markme[row, 'col']], rownames(df)[markme[row, 'row']],
          colnames(df)[markme[row, 'col']], rownames(df)[markme[row, 'row']])
    corrRect(p, namesMat = r) -> p
}
text(x=7.25, y=1.5, label=bquote(~log[10]~italic((p))), srt=270)
dev.off()
```

## Figure 4

### Correlation to other disorders

```r
library(ggplot2)
library(ggpubr)

corrs = readRDS('~/data/post_mortem_adhd/results/disorders_corrs_main.rds')
sources = unique(corrs$source)

# this leveling only affects the color ordering
dis_order = c('ASD', 'SCZ', 'BD', 'MDD', 'AAD', 'OCD')
col_labels = c('Autism Spectrum Disorder', 'Schizophrenia', 'Bipolar Disorder',
               'Major Depression Disorder', 'Alcohol Abuse or Dependence',
               'Obsessive Compulsive Disorder')
corrs$Disorders = factor(corrs$disorder,
                        levels=dis_order)

# just to share axis
ymax = .6
ymin = -.2
leg_size = 9
cap_size = 10
tick_size = 9
title_size = 11
# Bang Wong's color pallete
my_colors = c('#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2')

r = 'ACC' 
mycorrs = corrs[corrs$region == r, ]
mycorrs$id = sapply(1:nrow(mycorrs),
                  function(i) sprintf('%s [%d]',
                                      mycorrs[i, 'disorder'],
                                      which(sources == mycorrs[i, 'source'])))
# setting the order results appear in X axis
library(dplyr)
ranks = mycorrs %>% group_by(id) %>% summarize(Median = abs(median(corr, na.rm=TRUE)))
mylevels = c()
for (d in dis_order) {
    # assumes IDs start with the disorder
    myranks = ranks[grepl(ranks$id, pattern=paste0('^',d)),]
    mylevels = c(mylevels, myranks$id[order(myranks$Median, decreasing=T)])
}
# regroup levels based on the order we established before
mycorrs$xorder = factor(mycorrs$id, levels=mylevels)

p1 = ggplot(mycorrs, aes(x = xorder, y = corr, fill=Disorders)) +
    geom_violin(trim=FALSE) + 
    # fake continuous axis to add vertical lines later
    geom_line(aes(x = as.numeric(xorder), y=0), size = 1, color="#dc3220", alpha=0) + 
    # vertical lines separating disorders
    geom_vline(xintercept=c(4.5, 7.5, 10.5, 12.5),
               linetype="dashed", color = "grey", size=1) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=tick_size),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size=leg_size),
          legend.title = element_text(size=title_size),
          axis.text.y = element_text(size=tick_size),
          plot.title = element_text(size=title_size)) +
    ggtitle(r) + geom_hline(yintercept=0, linetype="dotted",
                                color = "#dc3220", size=1) +
   ylab(' ') + ylim(ymin, ymax) + 
   scale_fill_manual(breaks = levels(corrs$Disorder),
                     values = my_colors,
                     labels = col_labels,
                     drop=FALSE)

r = 'Caudate'
mycorrs = corrs[corrs$region == r, ]
mycorrs$id = sapply(1:nrow(mycorrs),
                  function(i) sprintf('%s [%d]',
                                      mycorrs[i, 'disorder'],
                                      which(sources == mycorrs[i, 'source'])))
# setting the order results appear in X axis
library(dplyr)
ranks = mycorrs %>% group_by(id) %>% summarize(Median = abs(median(corr, na.rm=TRUE)))
mylevels = c()
for (d in dis_order) {
    # assumes IDs start with the disorder
    myranks = ranks[grepl(ranks$id, pattern=paste0('^',d)),]
    mylevels = c(mylevels, myranks$id[order(myranks$Median, decreasing=T)])
}
# regroup levels based on the order we established before
mycorrs$xorder = factor(mycorrs$id, levels=mylevels)

p2 = ggplot(mycorrs, aes(x = xorder, y = corr, fill=Disorders)) +
    geom_violin(trim=FALSE) + 
    # fake continuous axis to add vertical lines later
    geom_line(aes(x = as.numeric(xorder), y=0), size = 1, color="#dc3220", alpha=0) + 
    # vertical lines separating disorders
    geom_vline(xintercept=c(2.5, 5.5, 8.5, 9.5, 10.5),
               linetype="dashed", color = "grey", size=1) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=tick_size),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size=leg_size),
          legend.title = element_text(size=title_size),
          axis.text.y = element_text(size=tick_size),
          plot.title = element_text(size=title_size)
          ) +
    ggtitle(r) + geom_hline(yintercept=0, linetype="dotted",
                                color = "#dc3220", size=1) +
   ylab(' ') + ylim(ymin, ymax) + 
   scale_fill_manual(breaks = levels(corrs$Disorder),
                     values = my_colors,
                     labels = col_labels,
                     drop=FALSE)

p = ggarrange(p1, p2, common.legend = T, legend='right', nrow=2, ncol=1,
          legend.grob=get_legend(p2)) 
p3 = annotate_figure(p, left = text_grob("Transcriptome correlation (rho)",
                                         rot = 90, size=cap_size))
```

### MAGMA results

```r
library(corrplot)
res = readRDS('~/data/post_mortem_adhd/results/MAGMA_main_res.rds')
res = res[res$POP == 'main', c('REGION', 'DISORDER', 'P')]
# we report ADHD in the text, and the other ones in the figure
res = res[res$DISORDER != 'ADHD', ]

plot_mat = matrix(nrow=2, ncol=length(unique(res$DISORDER)))
colnames(plot_mat) = unique(res$DISORDER)
rownames(plot_mat) = unique(res$REGION)
pvals = plot_mat
for (i in 1:nrow(res)) {
    plot_mat[res[i, 'REGION'], res[i, 'DISORDER']] = -log10(res[i, 'P'])
    pvals[res[i, 'REGION'], res[i, 'DISORDER']] = res[i, 'P']
}
quartz()
corrplot(plot_mat, is.corr=F, tl.col='black', p.mat = pvals,
         sig.level = .05, col.lim=c(0,8.5), cl.length=2,
         cl.align.text='l', cl.offset=.2, tl.cex = .8,
         col=colorRampPalette(c("white","#005AB5"))(200),
         mar = c(.5, .5, 1, 0))
text(x=9.15, y=1.5, label=bquote(~-log[10]~italic((p))), srt=270)
magma = recordPlot()

# adding some fake margin
library(cowplot)
fig4 = plot_grid(p3, magma, rel_heights = c(2.5, 1),
                 labels = c('A', 'B'), ncol=1) 
ggsave("~/data/post_mortem_adhd/figures/figure4.pdf", fig4, width=7.5, height=7.75,
       units="in")
```

# Supplemental figures

## SFigure 1

Diagram in Powerpoint, and added manually the figure from QoRTs.

## SFigure 2

Collage of 3 different views of ACC and Caudate, created using AFNI's SUMA, and
concatenated in Powerpoint.

## SFigure 3

### Raw PCA figures

```r
data = read.table('~/data/post_mortem_adhd/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = round(data)
sub_name = gsub(x=colnames(data), pattern='X', replacement='')
colnames(data) = sub_name
# this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
# its other replicate
data = data[, ! colnames(data) %in% c('66552')]

library(gdata)
df = read.xls('~/data/post_mortem_adhd/POST_MORTEM_META_DATA_JAN_2021.xlsx')
data = data[, colnames(data) %in% df$submitted_name]
df = df[df$submitted_name %in% colnames(data), ]
df = df[order(df$submitted_name), ]
data = data[, order(df$submitted_name)]

# cleaning up some variables
df$Individual = factor(df$hbcc_brain_id)
df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
df$MoD = factor(df$Manner.of.Death)
df$Sex = factor(df$Sex)
df$batch = factor(df$batch)
df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
df$SUB2 = 'no'
df[df$substance_group > 0, 'SUB2'] = 'yes'
df$SUB2 = factor(df$SUB2)
df$substance_group = factor(df$substance_group)
df$comorbid_group = factor(df$comorbid_group_update)
df$evidence_level = factor(df$evidence_level)
df$brainbank = factor(df$bainbank)
# replace the one subject missing population PCs by the median of their
# self-declared race and ethnicity
idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
med_pop = apply(df[idx, pop_pcs], 2, median)
df[which(is.na(df$C1)), pop_pcs] = med_pop
# combine brain bank and batch variables
df$BBB = factor(sapply(1:nrow(df),
                        function(x) sprintf('%s_%s',
                                    as.character(df[x,'brainbank']),
                                    as.character(df[x, 'batch']))))
df$BBB2 = NA                                                                        
df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
df[df$batch==3, 'BBB2'] = 2                                                         
df[df$batch==4, 'BBB2'] = 3      
df$BBB2 = factor(df$BBB2)
imWNH = which(df$C1 > 0 & df$C2 < -.075)
df$POP_BIN = 'other'
df[imWNH, 'POP_BIN'] = 'WNH'
df$POP_BIN = factor(df$POP_BIN)  
                    
library(GenomicFeatures)
txdb <- loadDb('~/one_drive/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
            "GENEID")
bt = read.csv('~/one_drive/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                tx_meta$TXCHROM != 'Y' &
                tx_meta$TXCHROM != 'MT')
data = data[imautosome, ]
tx_meta = tx_meta[imautosome, ]

library("DESeq2")
fm_str = '~ Diagnosis'
dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = df,
                                design = as.formula(fm_str))
vsd <- vst(dds, blind=FALSE) 
norm.cts <- assay(vsd) 
myvar = apply(norm.cts, 1, sd)

# similar to pcaExplorer, we select the top 300 genes in terms of variability
mytop = sort(myvar, index.return=T, decreasing = T)
top_data = norm.cts[mytop$ix[1:300],]
mypca <- prcomp(t(top_data), scale=TRUE) 
npcs = 12
pca_obj = mypca
var_explained_df <- data.frame(PC=factor(colnames(pca_obj$x)[1:npcs], 
                                         levels=colnames(pca_obj$x)[1:npcs]), 
                               var_explained=100*(pca_obj$sdev[1:npcs])^2/sum((pca_obj$sdev)^2))                             

# scree plot
library(ggplot2)
library(ggpubr)
p1 = ggplot(aes(x=PC, y=var_explained), data=var_explained_df)+ 
       geom_col()+ 
       geom_text(aes(label=round(var_explained, 1)),
                 position=position_dodge(width=0.9),
                 vjust=-0.5, size=2) + ylim(c(0, 80)) +
       geom_hline(yintercept=var_explained_df[7, 2], linetype="dashed",
                  color = "#dc3220") + 
       labs(title="Scree plot: RNAseq counts", 
            y='Percent variance explained (%)', 
            x='Principal components') +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))                                                                                        

# PC1 by PC2
plot_df = data.frame(PC1 = pca_obj$x[, 'PC1'], PC2 = pca_obj$x[, 'PC2'],
                     Region = df$Region)
plot_df$outlier = 'no'
idx = plot_df$PC1 > -10 & plot_df$PC1 < 10
plot_df[idx, ]$outlier = 'yes'
plot_df$outlier = as.factor(plot_df$outlier)
p2 = ggplot(aes(x=PC1, y=PC2, colour=outlier, fill=Region), data=plot_df)+ 
      geom_point(size = 3, shape = 21, stroke = .5) +
      labs(y = sprintf('PC2 (%.2f%%)', var_explained_df[2, 'var_explained']),
           x = sprintf('PC1 (%.2f%%)', var_explained_df[1, 'var_explained'])) +
      theme(legend.position="bottom") +
    scale_fill_manual(values=c('#FFC20A', "#0C7BDC")) +
    scale_color_manual(values=c('no'='black', 'yes'="red"), guide='none')
      
# PC2 by PC3
plot_df = data.frame(PC2 = pca_obj$x[, 'PC2'], PC3 = pca_obj$x[, 'PC3'],
                     BBB2 = df$BBB2)
plot_df$outlier = 'no'
# samples are in the same order as before
plot_df[idx, ]$outlier = 'yes'
plot_df$outlier = as.factor(plot_df$outlier)
p3 = ggplot(aes(x=PC2, y=PC3, fill=BBB2, colour=outlier), data=plot_df)+ 
      geom_point(size = 3, shape = 21, stroke = .5) +
      labs(y = sprintf('PC3 (%.2f%%)', var_explained_df[3, 'var_explained']),
           x = sprintf('PC2 (%.2f%%)', var_explained_df[2, 'var_explained']),
           color = 'Brain bank / batch') +
      theme(legend.position="bottom") +
      scale_fill_manual(values=c('#0C7BDC', "#F0E442", "#009E73")) + 
      scale_color_manual(values=c('no'='black', 'yes'="red"), guide='none')
```

### Variance partition plot

Let's focus on the clean data here, so we can show better the contribution of
the different covariates:

```r
data = read.table('~/data/post_mortem_adhd/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = round(data)
sub_name = gsub(x=colnames(data), pattern='X', replacement='')
colnames(data) = sub_name
# this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
# its other replicate
data = data[, ! colnames(data) %in% c('66552')]
# outliers based on PCA plots
outliers = c('68080','68096', '68108', '68084', '68082')
data = data[, ! colnames(data) %in% outliers]

library(gdata)
df = read.xls('~/data/post_mortem_adhd/POST_MORTEM_META_DATA_JAN_2021.xlsx')
# sync data and df
data = data[, colnames(data) %in% df$submitted_name]
df = df[df$submitted_name %in% colnames(data), ]
df = df[order(df$submitted_name), ]
data = data[, order(df$submitted_name)]

# cleaning up some variables
df$Individual = factor(df$hbcc_brain_id)
df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
df$MoD = factor(df$Manner.of.Death)
df$Sex = factor(df$Sex)
df$batch = factor(df$batch)
df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
df$SUB2 = 'no'
df[df$substance_group > 0, 'SUB2'] = 'yes'
df$SUB2 = factor(df$SUB2)
df$substance_group = factor(df$substance_group)
df$comorbid_group = factor(df$comorbid_group_update)
df$evidence_level = factor(df$evidence_level)
df$brainbank = factor(df$bainbank)
# replace the one subject missing population PCs by the median of their
# self-declared race and ethnicity
idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
med_pop = apply(df[idx, pop_pcs], 2, median)
df[which(is.na(df$C1)), pop_pcs] = med_pop
# combine brain bank and batch variables
df$BBB = factor(sapply(1:nrow(df),
                        function(x) sprintf('%s_%s',
                                    as.character(df[x,'brainbank']),
                                    as.character(df[x, 'batch']))))
df$BBB2 = NA                                                                        
df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
df[df$batch==3, 'BBB2'] = 2                                                         
df[df$batch==4, 'BBB2'] = 3      
df$BBB2 = factor(df$BBB2)
imWNH = which(df$C1 > 0 & df$C2 < -.075)
df$POP_BIN = 'other'
df[imWNH, 'POP_BIN'] = 'WNH'
df$POP_BIN = factor(df$POP_BIN)  

# removing non-autosome genes
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem_adhd/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
            "GENEID")
bt = read.csv('~/data/post_mortem_adhd/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                tx_meta$TXCHROM != 'Y' &
                tx_meta$TXCHROM != 'MT')
data = data[imautosome, ]
tx_meta = tx_meta[imautosome, ]

# remove constant genes (including zeros) as it breaks PCA
const_genes = apply(data, 1, sd) == 0
data = data[!const_genes, ]

library("DESeq2")
# making sure any numeric covariates are scaled
num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
        'C1', 'C2', 'C3', 'C4', 'C5')
for (var in num_vars) {
    df[, var] = scale(df[, var])
}

# creating DESeq2 object
fm_str = '~ RINe + C1 + BBB2 + comorbid_group + SUB2 + Diagnosis'
dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = df,
                                design = as.formula(fm_str))
# remove genes based on how many subjects have zero counts
# the minimum number of subjects with zero counts in a gene for the gene
# to be kept corresponds to the smallest group size (number of Cases)
min_subjs = min(table(df$Diagnosis))
keep <- rowSums(counts(dds) == 0) <= min_subjs
dds <- dds[keep,]

# remove genes based on filterByExpr()
library(edgeR)
design = model.matrix(as.formula(fm_str), data=colData(dds))
isexpr <- filterByExpr(counts(dds), design=design)
ddsExpr = dds[isexpr, ]

# variancePartition plot
library(variancePartition)
vsd <- vst(ddsExpr, blind=FALSE) 
norm.cts <- assay(vsd) 
myvar = apply(norm.cts, 1, sd)

form <- ~ (C1 + C2 + C3 + C4 + C5 + Age + RINe + (1|BBB2) + (1|comorbid_group) +
           (1|SUB2) + (1|Diagnosis) + PMI + (1|MoD) + (1|Sex) +
           (1|evidence_level) + (1|Region))
# remove genes with zero variance
norm.cts.var = norm.cts[myvar>0,]
# varPart <- fitExtractVarPartModel(norm.cts.var, form, df)
# ran it in the cluster
load('~/tmp/varpart.rdata')
colnames(varPart) = c("Brain bank / batch", "Comorbidities", "Diagnosis", 
                      "Evidence level", "Mode of death", "Brain region", "Sex",
                      "Substance abuse", "C1", "C2", "C3", "C4", "C5", "Age",
                      "RINe", "Post-mortem interval", "Residuals")
vp <- sortCols( varPart )
fig = plotVarPart( vp )
```

Combining all plots:

```r
fig_mod = fig + theme(axis.text.x = element_text(size = 9, angle = 55),
                      axis.title.y = element_text(size = 9),
                      axis.text.y = element_text(size = 8),
                      plot.margin = margin(l=7, b=0, r=7, t=15))
p = ggarrange(fig_mod, p1, p2, p3, common.legend=F, ncol=2, nrow=2,
              labels='AUTO')

ggsave("~/data/post_mortem_adhd/figures/sfigure3.pdf", p, width=6.5, height=6,
       units="in")
```

## SFigure 4

```r
get_data_and_df = function(myregion = NA) {
    data = read.table('~/data/post_mortem_adhd/adhd_rnaseq_counts.txt', header=1)
    rownames(data) = data[,1]
    data[,1] = NULL
    data = round(data)
    sub_name = gsub(x=colnames(data), pattern='X', replacement='')
    colnames(data) = sub_name
    # removing same outliers we did for DGE analysis
    data = data[, ! colnames(data) %in% c('66552')]
    # outliers based on PCA plots
    outliers = c('68080','68096', '68108', '68084', '68082')
    data = data[, ! colnames(data) %in% outliers]

    library(gdata)
    df = read.xls('~/data/post_mortem_adhd/POST_MORTEM_META_DATA_JAN_2021.xlsx')
    data = data[, colnames(data) %in% df$submitted_name]
    df = df[df$submitted_name %in% colnames(data), ]
    df = df[order(df$submitted_name), ]
    data = data[, order(df$submitted_name)]

    if (is.na(myregion)) {
        # only a single entry per individual
        keep_me = !duplicated(df$hbcc_brain_id)
    } else {
        keep_me = df$Region == myregion
    }
    data = data[, keep_me]
    df = df[keep_me, ]

    # cleaning up some variables
    df$Individual = factor(df$hbcc_brain_id)
    df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
    df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
    df$MoD = factor(df$Manner.of.Death)
    df$Sex = factor(df$Sex)
    df$batch = factor(df$batch)
    df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
    df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
    df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
    df$SUB2 = 'no'
    df[df$substance_group > 0, 'SUB2'] = 'yes'
    df$SUB2 = factor(df$SUB2)
    df$substance_group = factor(df$substance_group)
    df$comorbid_group = factor(df$comorbid_group_update)
    df$evidence_level = factor(df$evidence_level)
    df$brainbank = factor(df$bainbank)
    # replace the one subject missing population PCs by the median of their
    # self-declared race and ethnicity
    idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
    pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
    med_pop = apply(df[idx, pop_pcs], 2, median)
    df[which(is.na(df$C1)), pop_pcs] = med_pop
    df$BBB = factor(sapply(1:nrow(df),
                            function(x) sprintf('%s_%s',
                                        as.character(df[x,'brainbank']),
                                        as.character(df[x, 'batch']))))
    df$BBB2 = NA                                                                        
    df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
    df[df$batch==3, 'BBB2'] = 2                                                         
    df[df$batch==4, 'BBB2'] = 3      
    df$BBB2 = factor(df$BBB2)
    imWNH = which(df$C1 > 0 & df$C2 < -.075)
    df$POP_BIN = 'other'
    df[imWNH, 'POP_BIN'] = 'WNH'
    df$POP_BIN = factor(df$POP_BIN)        
    # df$RINc = cut(df$RINe, breaks = 4)  
    # bining so DESeq2 can do its own filyering automatically
    breaks = quantile(df$RINe, probs = seq(0, 1, by = 0.25))
    df$RINc = cut(df$RINe, breaks=breaks, labels=c('q1', 'q2', 'q3', 'q4'),
                include.lowest=T)

    library(GenomicFeatures)
    txdb <- loadDb('~/data/post_mortem_adhd/Homo_sapies.GRCh38.97.sqlite')
    txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
                "GENEID")
    bt = read.csv('~/data/post_mortem_adhd/Homo_sapiens.GRCh38.97_biotypes.csv')
    bt_slim = bt[, c('gene_id', 'gene_biotype')]
    bt_slim = bt_slim[!duplicated(bt_slim),]
    txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
    tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
    tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
    imautosome = which(tx_meta$TXCHROM != 'X' &
                    tx_meta$TXCHROM != 'Y' &
                    tx_meta$TXCHROM != 'MT')
    data = data[imautosome, ]
    tx_meta = tx_meta[imautosome, ]

    # remove constant genes (including zeros) as it breaks PCA
    const_genes = apply(data, 1, sd) == 0
    data = data[!const_genes, ]

    min_subjs = min(table(df$Diagnosis))
    keep <- rowSums(data == 0) <= min_subjs
    data <- data[keep,]

    res = list(df=df, data=data)
    return(res)
}

m = get_data_and_df(NA)

# run nonparametric t-tests for numeric variables
num_vars = c('Age', 'PMI', 'C1', 'C2', 'C3', 'C4', 'C5', 'RINe')
mypvals = c()
mystats = c()
for (x in num_vars) {
    res = wilcox.test(as.formula(sprintf('%s ~ Diagnosis', x)), data=m$df)
    mypvals = c(mypvals, res$p.value)
    mystats = c(mystats, res$statistic)
}

categ_vars = c('MoD', 'SUB2', 'comorbid_group', 'Sex', 'evidence_level', 'BBB2')
for (x in categ_vars) {
    res = chisq.test(table(m$df$Diagnosis, m$df[, x]))
    mypvals = c(mypvals, res$p.value)
    mystats = c(mystats, res$statistic)
}

myvars = c(num_vars, categ_vars)
DX_pvals = mypvals
DX_plot = -log10(mypvals) * sign(mystats)
names(DX_plot) = myvars
```

Now we run it for the PCA of ACC and Caudate separately. Make sure we save the
data for the scree plot as well:

```r
get_correlations = function(data, df) {
    # checking which PCs are associated with our potential nuiscance variables
    set.seed(42)
    mypca <- prcomp(t(data), scale=TRUE)
    # how many PCs to keep... using Kaiser thredhold, close to eigenvalues < 1
    library(nFactors)
    eigs <- mypca$sdev^2
    nS = nScree(x=eigs)
    keep_me = seq(1, nS$Components$nkaiser)
    mydata = data.frame(mypca$x[, keep_me])
    # create main metadata data frame including metadata and PCs
    data.pm = cbind(df, mydata)
    rownames(data.pm) = df$hbcc_brain_id
    cat('Using', nS$Components$nkaiser, 'PCs from possible', ncol(data), '\n')

    # check which PCs are associated at nominal p<.01
    pc_vars = colnames(mydata)
    num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                        dimnames=list(num_vars, pc_vars))
    num_pvals = num_corrs
    for (x in num_vars) {
        for (y in pc_vars) {
            res = cor.test(data.pm[, x], data.pm[, y], method='spearman')
            num_corrs[x, y] = res$estimate
            num_pvals[x, y] = res$p.value
        }
    }

    categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                            dimnames=list(categ_vars, pc_vars))
    categ_pvals = categ_corrs
    for (x in categ_vars) {
        for (y in pc_vars) {
            res = kruskal.test(data.pm[, y], data.pm[, x])
            categ_corrs[x, y] = res$statistic
            categ_pvals[x, y] = res$p.value
        }
    }
    mypvals = rbind(categ_pvals, num_pvals)
    mycorrs = rbind(categ_corrs, num_corrs)

    res = list(pvals=mypvals, corrs=mycorrs, pca=mypca)
    return(res)
}

m = get_data_and_df('ACC')
ACC_res = get_correlations(m$data, m$df)
m = get_data_and_df('Caudate')
Caudate_res = get_correlations(m$data, m$df)
```

Now, let's work on plotting them. Scree first:

```r
library(ggplot2)

pca_obj = ACC_res$pca
npcs = 12
var_explained_df <- data.frame(PC=factor(colnames(pca_obj$x)[1:npcs],
                                         levels=colnames(pca_obj$x)[1:npcs]),
                               var_explained=100*(pca_obj$sdev[1:npcs])^2/sum((pca_obj$sdev)^2))
p1 = ggplot(aes(x=PC, y=var_explained), data=var_explained_df)+
  geom_col()+
  geom_hline(yintercept=var_explained_df[7, 2], linetype="dashed", color = "#dc3220") +
  labs(title="ACC",
        y='Percent variance explained (%)',
        x='Principal components') +
   ylim(0, 50) + geom_text(aes(label=round(var_explained, 1)),
                 position=position_dodge(width=0.9),
                 vjust=-0.5, size=2) +
   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))

pca_obj = Caudate_res$pca
var_explained_df <- data.frame(PC=factor(colnames(pca_obj$x)[1:npcs],
                                         levels=colnames(pca_obj$x)[1:npcs]),
                               var_explained=100*(pca_obj$sdev[1:npcs])^2/sum((pca_obj$sdev)^2))
p2 = ggplot(aes(x=PC, y=var_explained), data=var_explained_df)+
  geom_col()+
  geom_hline(yintercept=var_explained_df[8, 2], linetype="dashed", color = "#dc3220") +
  labs(title="Caudate",
        y=' ',
        x=' ') +
   ylim(0, 50) + geom_text(aes(label=round(var_explained, 1)),
                 position=position_dodge(width=0.9),
                 vjust=-0.5, size=2) +
   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
```

Now the correlation plot:

```r
library(corrplot)
ACC_plot = -log10(ACC_res$pvals) * sign(ACC_res$corrs)
cnames = sapply(colnames(ACC_plot), function(x) sprintf('ACC: %s', x))
plot_df = t(ACC_plot)
rownames(plot_df) = cnames
Caudate_plot = -log10(Caudate_res$pvals) * sign(Caudate_res$corrs)
cnames = sapply(colnames(Caudate_plot), function(x) sprintf('Caudate: %s', x))
colnames(Caudate_plot) = cnames
plot_df = rbind(plot_df, t(Caudate_plot))
junk = t(data.frame(DX_plot))
plot_df2 = rbind(plot_df, junk[, match(colnames(plot_df), colnames(junk))])
rownames(plot_df2)[nrow(plot_df2)] = 'Diagnosis'
mylim = max(abs(plot_df2))

colnames(plot_df2) = c("Mode of death", "Substance abuse", "Comorbidities",
                       "Sex", "Evidence level", "Brain bank / batch", "Age",
                       "Post-mortem interval", "C1", "C2", "C3", "C4", "C5",
                       "RINe")

# select these variables in the plot
mypvals = ACC_res$pvals
rownames(mypvals) = colnames(plot_df2)
ncomps = ncol(mypvals) * nrow(mypvals)
print(unique(rownames(which(mypvals <= .05/ncomps, arr.ind = T))))
mypvals = Caudate_res$pvals
rownames(mypvals) = colnames(plot_df2)
ncomps = ncol(mypvals) * nrow(mypvals)
print(unique(rownames(which(mypvals <= .05/ncomps, arr.ind = T))))
ncomps = length(DX_pvals)
print(names(DX_plot)[which(DX_pvals <= .05/ncomps, arr.ind = T)])

quartz()
par(xpd=TRUE)
corrplot(plot_df2, is.corr=F, col.lim=c(-mylim, mylim), tl.col='black',
         tl.cex=.8, cl.cex=.8, cl.align.text = 'l', cl.offset=.25,
         mar=c(1,0,3.5,0),
         col=colorRampPalette(c("#dc3220", "white","#005AB5"))(200))
text(x=18.5, y=9, label=bquote(~-log[10]~italic((p))), srt=-270)
width = 1
ybottom = 16.65
xpos = c(1.5, 2.5, 5.5, 8.5, 13.5)
ytop = c(6.3, 5.1, 6.6, 1.1, 2)
rect(xpos, ybottom, xpos+width, ybottom+ytop, border='#dc3220', lty='dashed', lwd=1.5)
corrs = recordPlot()
```

Putting it all together:

```r
library(cowplot)
bottom = plot_grid(p1, p2, rel_widths = c(1, 1), ncol=2) 
sfig4 = plot_grid(corrs, bottom, rel_heights = c(1.5,1), ncol=1, labels='AUTO')
ggsave("~/data/post_mortem_adhd/figures/sfigure4.pdf", sfig4, width=7.5, height=7.75,
       units="in")
```


## SFigure 5

```r
library(gdata)
library(ggplot2)

df = read.xls('~/data/post_mortem_adhd/POST_MORTEM_META_DATA_JAN_2021.xlsx')
# removing the same samples I did for DGE
df = df[df$submitted_name != '66552',]
outliers = c('68080','68096', '68108', '68084', '68082')
df = df[! df$submitted_name %in% outliers,]

df = df[!duplicated(df$hbcc_brain_id),]
df$POP_CODE = 'White non-Hispanic'
idx = which(df$Ethnicity.x == 'Hispanic' & df$Race.x == 'White')
df[idx, 'POP_CODE'] = 'White Hispanic'
idx = which(df$Race.x == 'African American' | df$Race.x == 'Black or African-American')
df[idx, 'POP_CODE'] = 'African-American'

df$POP_CODE = factor(df$POP_CODE)

p = ggplot(df, aes(x=C1, y=C2, color=POP_CODE)) + geom_point() +
    geom_vline(aes(xintercept=0), color='#dc3220', linetype='dashed') +
    geom_hline(aes(yintercept=-.075), color='#dc3220', linetype='dashed') +
    labs(color='Sub-populations') +
    scale_colour_manual(values=c('#0C7BDC', "#F0E442", "#009E73"))

ggsave("~/data/post_mortem_adhd/figures/sfigure5.pdf", p, width=5, height=5,
        units="in")
```

## SFigure 6

```r
load('~/data/post_mortem_adhd/results/main_DGE.RData')
library(DESeq2)

get_plot_list = function(dds) {
    res = results(dds, name = "Diagnosis_Case_vs_Control", alpha=.05)
    res = res[order(res$pvalue), ]
    gene_list = rownames(res)[which(round(res$padj, 2) <= .05)]

    trans.dds <- vst(dds, blind=FALSE)
    norm.cts <- assay(trans.dds)

    covars = model.matrix(~ RINe + C1 + BBB2 + comorbid_group + SUB2,
                        data=colData(dds))
    dsn = model.matrix(~ Diagnosis, data=colData(dds))
    mat <- limma::removeBatchEffect(norm.cts, covariates=covars, design=dsn)

    gnames = data.frame(full=rownames(counts(dds)),
                        nov=substr(rownames(counts(dds)), 1, 15))
    mart = readRDS('~/data/post_mortem_adhd/mart_rnaseq.rds')
    gnames = merge(gnames, mart, by.x='nov', by.y='ensembl_gene_id')
    keep_me = gnames$full %in% gene_list
    gene_ids = gnames[keep_me, ]

    resid_expr = reshape2::melt(mat[gene_ids$full,])
    if (nrow(gene_ids) > 1) {
        colnames(resid_expr) = c('gene', 'submitted_name', 'normCount')
    } else {
        resid_expr$gene = gene_ids$full
        resid_expr$submitted_name = rownames(resid_expr)
        resid_expr$normCount = resid_expr$value
    }
    junk = colData(trans.dds)[, c('Diagnosis', 'submitted_name')]
    resid_expr = merge(resid_expr, junk, by='submitted_name')
    resid_expr = merge(resid_expr, gene_ids, by.x='gene', by.y='full')

    # plotting each of the significant genes
    library(ggpubr)
    library(ggbeeswarm)
    myplots = list()
    clrs = c("#1AFF1A", "#4B0092")
    for (g in 1:nrow(gene_ids)) {
        hgnc = gene_ids[g, 'hgnc_symbol']
        if (hgnc == '') {
            hgnc = gene_ids[g, 'nov']
        }
        cat(gene_ids[g, 'nov'], hgnc, '\n')
        d = as.data.frame(resid_expr[resid_expr$gene == gene_list[g],])
        p = (ggplot(d, aes(x=Diagnosis, y=normCount, color = Diagnosis,
                            fill = Diagnosis)) + 
                scale_y_log10() +
                geom_boxplot(alpha = 0.4, outlier.shape = '+', width = 0.8,
                            lwd = 0.5, notch=F) +
                scale_color_manual(values = clrs, labels=c('Unaffected', 'ADHD')) +
                scale_fill_manual(values = clrs, labels=c('Unaffected', 'ADHD')) +
                theme_bw() + theme(axis.text.x = element_blank(),
                                axis.title.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                axis.title.y = element_blank(),
                                axis.text.y = element_text(size = 10),
                                plot.title = element_text(size = 11)) +
                ggtitle(hgnc))
        myplots[[g]] = p + theme(legend.position = "none")
    }
    return(myplots)
}

library(ggpubr)
p.acc = ggarrange(plotlist=get_plot_list(dds.ACC), common.legend=T,
                  legend='none')
p.caudate = get_plot_list(dds.Caudate)[[1]] +
            ylab('Normalized\ngene counts') +
            theme(legend.position = "right",
                  axis.title.y = element_text(size=11, angle=90))

library(cowplot)
bottom = plot_grid(p.caudate, NULL, rel_widths = c(1, 1.25), ncol=2) 
sfig6 = plot_grid(NULL, p.acc, NULL, bottom, rel_heights = c(.15,3.75,.15,1), ncol=1,
                  labels=c('A. Anterior cingulate cortex', '', 'B. Caudate', ''),
                  label_x = 0, hjust=0)
ggsave("~/data/post_mortem_adhd/figures/sfigure6.pdf", sfig6, width=7.5, height=9,
       units="in")
```

## SFigure 7

```r
library(ggpubr)
tsize = 12
ysize = 10
xsize = 9
msize = 3

df = read.csv('~/data/post_mortem_adhd/results/GSEA_main_Caudate_geneontology_Cellular_Component_noRedundant_10K.csv')
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore', 'FDR')]
df = df[round(df$FDR, 2) <= .05, ]
df[df$FDR == 0, 'FDR'] = 1e-4

df$description = factor(df$description,
                        levels=df$description[sort(df$FDR,
                                                   index.return=T,
                                                   decreasing=T)$ix])
color_me = c('neuron spine', 'synaptic membrane', 'glutamatergic synapse',
             'postsynaptic specialization', 'neuron to neuron synapse',
             'Schaffer collateral - CA1 synapse', 'axon part',
             'GABA-ergic synapse', 'presynapse',
             'neuron projection terminus', 'dendritic shaft',
             'neuronal cell body', 'neuromuscular junction', 'coated membrane',
             'myelin sheath', 'excitatory synapse')
label_colors <- ifelse(levels(df$description) %in% color_me, "#dc3220", "grey20")
df$Behavior = ifelse(df$description %in% color_me, "neuro", "notneuro")

p <- ggplot(df, aes(y=-log10(FDR), x=description, fill=Behavior)) +
  geom_dotplot(dotsize=msize, binaxis='y', stackdir='center') + coord_flip() +
  geom_hline(yintercept=-log10(.1), linetype="dashed",
                                color = "black", size=1) + theme(legend.position="bottom") +
    geom_hline(yintercept=-log10(.05), linetype="dotted",
                                color = "black", size=1) + theme(legend.position="bottom") +
    theme(axis.text.y = element_text(colour = label_colors, size=ysize),
          axis.title.y = element_blank(),
          plot.title = element_text(size = tsize),
          axis.text.x = element_text(size=xsize),
          axis.title.x = element_blank(),
          plot.margin = margin(l=3, b=4, r=3, t=7)) +
    ylim(0, 5) +
    scale_fill_manual(values=c("#DC3220", 'grey60'))

p.cau = p + ggtitle('Caudate') + theme(axis.title.x = element_text(size=xsize)) +
        ylab(bquote(~-log[10]~italic((p[adjusted]))))

df = read.csv('~/data/post_mortem_adhd/results/GSEA_main_ACC_geneontology_Cellular_Component_noRedundant_10K.csv')
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore', 'FDR')]
df = df[round(df$FDR, 2) <= .05, ]
df[df$FDR == 0, 'FDR'] = 1e-4

df$description = factor(df$description,
                        levels=df$description[sort(df$FDR,
                                                   index.return=T,
                                                   decreasing=T)$ix])

color_me = c('nothing')
label_colors <- ifelse(levels(df$description) %in% color_me, "#dc3220", "grey20")
df$Behavior = ifelse(df$description %in% color_me, "neuro", "notneuro")
df$Behavior = factor(df$Behavior, levels=c('neuro', 'notneuro'))

p <- ggplot(df, aes(y=-log10(FDR), x=description, fill=Behavior, size=msize)) +
  geom_dotplot(dotsize=msize*1.2, binaxis='y', stackdir='center') + coord_flip() +
  geom_hline(yintercept=-log10(.1), linetype="dashed",
                                color = "black", size=1) + theme(legend.position="bottom") +
    geom_hline(yintercept=-log10(.05), linetype="dotted",
                                color = "black", size=1) + theme(legend.position="bottom") +
    theme(axis.text.y = element_text(colour = label_colors, size=ysize),
          axis.title.y = element_blank(),
          plot.title = element_text(size = tsize),
          axis.text.x = element_text(size=xsize),
          axis.title.x = element_blank(),
          plot.margin = margin(l=0, b=360, r=10, t=7)) +
    scale_fill_manual(values=c("#DC3220", 'grey60'), drop=F) + ylim(0, 5)

p.acc = p + ggtitle('ACC') + theme(axis.title.x = element_text(size=xsize)) +
        ylab(bquote(~-log[10]~italic((p[adjusted]))))

p2 = ggarrange(p.cau, p.acc, common.legend=T, ncol=2, nrow=1, widths=c(1.2,1),
               legend='none', labels='AUTO')
ggsave("~/data/post_mortem_adhd/figures/sfigure7.pdf", p2, width=7.5, height=8,
        units="in")
```

## SFigure 8

The semantic space plot is done in R, using a modified version of the code we
got from REVIGO. The other two plots are screenshots from the REVIGO website
results. We get these results after pasting the 40 significant FDR q < .05
Cellular Component results in the Caudate. I also selected the resulting list ot
be large (0.9), that the values represent P-values, to remove obsolete GO
terms, working with the Homo Sapiens database, and using the Lin similarity
metric. At the time, the Gene Ontology database was Tuesday, March 22, 2022.

```r
library( ggplot2 );
library( scales );

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","value","uniqueness","dispensability");
# removing obsolete tag as it was not present in the sets I used in WebGestalt
revigo.data <- rbind(c("GO:0099572","postsynaptic specialization",1.80789765819105,-3.88400328988402,-5.13739105957345,2.53529412004277,-300,0.736969043926848,0),
c("GO:1990351","transporter complex",2.19379394195697,4.08155640787265,-4.18246512368601,2.61909333062674,-3.36110658801899,0.816474409214203,0),
c("GO:0043209","myelin sheath",0.248453771739705,-4.81815464422998,3.88738081691802,1.68124123737559,-1.75848170961614,0.997739294793665,0.00162244),
c("GO:0070469","respirasome",0.544483797642332,5.48839850266622,4.46433184821522,2.01703333929878,-2.31913395602603,0.93714644088994,0.00176041),
c("GO:0030427","site of polarized growth",0.914521330020616,-6.84683215766251,-0.312015193100905,2.2405492482826,-2.34756740692001,0.997439805098722,0.00186525),
c("GO:0044309","neuron spine",0.930380081408257,-0.422099585751614,5.76944165649107,2.24797326636181,-4.66154350639539,0.918704714766417,0.00186894),
c("GO:0031252","cell leading edge",2.2307976951948,-5.47753348665079,0.786730816469546,2.62634036737504,-1.31922089869559,0.997186490647854,0.00207805),
c("GO:0036019","endolysosome",0.153301263413861,-5.90180482717856,2.32009615359036,1.47712125471966,-1.32604891992159,0.974123199620079,0.05689351),
c("GO:0048475","coated membrane",0.49690754347941,3.86349509856744,5.14403040621447,1.97772360528885,-1.88777856001282,0.93763517257539,0.1263669),
c("GO:0019898","extrinsic component of membrane",1.68102764708992,3.42697431584769,3.80976958233819,2.50379068305718,-1.38086235480548,0.930552428983581,0.14292921),
c("GO:0070069","cytochrome complex",0.22202251942697,3.88854636059026,-6.24145906942097,1.63346845557959,-1.68780937365379,0.849154291190375,0.23867133),
c("GO:0016469","proton-transporting two-sector ATPase complex",0.280171274514987,4.93450992785269,-1.04416832187025,1.73239375982297,-2.35571181473037,0.773324243333528,0.24439681),
c("GO:1903293","phosphatase complex",0.311888777290268,6.04498773855961,-3.88439873844692,1.77815125038364,-2.66395013729136,0.787165170400544,0.24712996),
c("GO:0030684","preribosome",0.401755035153566,4.92208514900057,-5.56146900929428,1.88649072517248,-1.8887826324559,0.821085472203846,0.25383178),
c("GO:0120114","Sm-like protein family complex",0.449331289316488,3.69427506836137,-5.26615438677703,1.93449845124357,-1.55317541765007,0.840277590874042,0.25691137),
c("GO:0098798","mitochondrial protein-containing complex",1.58058888830153,6.20237585629081,-1.49946311568005,2.47712125471966,-2.81703868616152,0.778277465779588,0.29747178),
c("GO:1905360","GTPase complex",0.195591267114236,5.83434651329172,-4.44668124063308,1.57978359661681,-1.42042166745485,0.794505608925383,0.4192086),
c("GO:0042611","MHC protein complex",0.13744251202622,3.52614855686166,-1.52598070797693,1.43136376415899,-1.76911651623963,0.776933544716485,0.42126077),
c("GO:0031201","SNARE complex",0.264312523127346,4.53453445921042,-1.23862990331234,1.70757017609794,-1.36967306185932,0.770932843236051,0.44457967),
c("GO:0030964","NADH dehydrogenase complex",0.269598773589893,4.84675860273709,-2.2621298870761,1.7160033436348,-1.38638522237082,0.719251499643988,0.44532617),
c("GO:1990204","oxidoreductase complex",0.681926309668552,5.21305950297771,-3.71555939489758,2.11394335230684,-1.96815874976306,0.773795736366881,0.46787318),
c("GO:0031514","motile cilium",1.30041761378654,-0.77245884020487,5.84466177589041,2.39269695325967,-2.00759102068833,0.907302520762106,0.47072525),
c("GO:0043198","dendritic shaft",0.195591267114236,-1.34775038182477,5.134630977349,1.57978359661681,-2.53526356806634,0.89419268938406,0.47986935),
c("GO:1904949","ATPase complex",0.507480044404504,5.66548223817041,-3.60245883782271,1.98677173426624,-1.55571520554076,0.779006746679486,0.49004777),
c("GO:1905368","peptidase complex",0.639636305968177,5.28826132212523,-4.12790553186024,2.08635983067475,-1.51559535881801,0.774941757548856,0.50134433),
c("GO:0060076","excitatory synapse",0.27488502405244,-3.31862529769077,-6.09286709491123,1.72427586960079,-1.31851534546173,0.790455755302001,0.53310826),
c("GO:0044306","neuron projection terminus",0.692498810593646,-0.135890411702787,5.87721100190712,2.12057393120585,-2.55991433991181,0.920200063950927,0.54273899),
c("GO:0031594","neuromuscular junction",0.370037532378284,-3.66482116689501,-5.98010111487439,1.85125834871908,-2.20881845436474,0.785198438025596,0.5495941),
c("GO:0098685","Schaffer collateral - CA1 synapse",0.380610033303378,-4.18069751182532,-4.76804726141264,1.86332286012046,-3.78707498780825,0.784687811879604,0.55120953),
c("GO:0098982","GABA-ergic synapse",0.385896283765925,-4.29408196216163,-5.12065847646183,1.86923171973098,-2.70789440672693,0.784436997242679,0.55200396),
c("GO:0005743","mitochondrial inner membrane",2.62198022942327,5.31927771013672,2.9521275122047,2.69635638873333,-1.3267316984797,0.877328170975456,0.62435707),
c("GO:0043025","neuronal cell body",2.61140772849818,-2.4050586095984,3.14301776659146,2.69460519893357,-2.23032645023267,0.926780204649808,0.62886877),
c("GO:0098978","glutamatergic synapse",1.71803140032775,-3.58406770454158,-5.59968229872664,2.51321760006794,-300,0.753971074718595,0.65406465),
c("GO:0005681","spliceosomal complex",1.0361050906592,6.4508626626043,-2.78771848600848,2.29446622616159,-1.56195779701283,0.789963859725258,0.65424163),
c("GO:0098984","neuron to neuron synapse",1.85547391235397,-3.96848993542213,-5.5369867951831,2.54654266347813,-300,0.752216076885374,0.66035681),
c("GO:0097060","synaptic membrane",1.99291642438019,-2.69636461111216,-3.88740061518261,2.57749179983723,-4.5654310959658,0.683478320852503,0.66849862),
c("GO:0044298","cell body membrane",0.169160014801501,-1.28951321411989,1.73471228992195,1.51851393987789,-1.54915820710792,0.897746016882989,0.70169915),
c("GO:0098793","presynapse",2.7488502405244,-3.72677685005849,-5.35769416488586,2.71683772329952,-2.66870068229318,0.74297970694927,0.70348715));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
# data$value is already log10(pvalue)!
one.data$value <- -as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );

#p-value for 0 shouldn't be -300
one.data$value[one.data$value == 300] = -log10(1e-4)
res = kmeans(one.data[, c('plot_X', 'plot_Y')], centers=4, algorithm='Forgy') 
one.data$cluster = factor(res$cluster)
p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour=cluster, size = plot_size), alpha = I(0.6), show.legend=F) + scale_size_area();
p1 <- p1 + scale_colour_manual(values=c('#E69F00', '#56B4E9', '#009E73', '#F0E442'))
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); 
ex <- one.data [ one.data$dispensability < 0.15, ];
# force legend a bit to the right
ex[ex$description == 'site of polarized growth', 'plot_X'] = -5.5
p1 <- p1 + geom_label( data = ex, aes(plot_X, plot_Y, label = description),
                      colour = I(alpha("black", 1)), size = 3 , fill='white');
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank(),
                 axis.text.x = element_text(size=10),
                 axis.text.y = element_text(size=10),
                 axis.title.x = element_text(size=11),
                 axis.title.y = element_text(size=11)) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
p1 = p1 + labs(size = bquote(~-log[10]~italic(p[adjusted])))

ggsave("~/tmp/semantic.pdf", p1, width=6.5, height=4.5, units="in")
```

## SFigure 9

```r
library(corrplot)

make_corrplot = function(g) {
    r = 'ACC'
    file_name = sprintf('~/data/post_mortem_adhd/results/GSEA_%s_%s_developmental_10K.csv',
                        g, r)
    res = read.csv(file_name)
    res = res[order(res$geneSet), c('link', 'normalizedEnrichmentScore', 'pValue')]
    dev = res
    dev$Region = r

    df = matrix(nrow = 2, ncol = 6, dimnames=list(c('ACC', 'Caudate'),
                                                c('pandev', res$link[1:5])))
    pvals = df
    i = 1
    for (j in 1:ncol(df)) {
        idx = dev$Region == rownames(df)[i] & dev$link == colnames(df)[j]
        if (dev[idx, 'pValue'] == 0) {
            dev[idx, 'pValue'] = 1e-4
        }
        df[i, j] = (sign(dev[idx, 'normalizedEnrichmentScore']) *
                    -log10(dev[idx, 'pValue']))
        pvals[i, j] = dev[idx, 'pValue']
    }

    r = 'Caudate'
    file_name = sprintf('~/data/post_mortem_adhd/results/GSEA_%s_%s_developmental_10K.csv',
                        g, r)
    res = read.csv(file_name)
    res = res[order(res$geneSet), c('link', 'normalizedEnrichmentScore', 'pValue')]
    res$Region = r
    dev = res

    i = 2
    for (j in 1:ncol(df)) {
        idx = dev$Region == rownames(df)[i] & dev$link == colnames(df)[j]
        if (dev[idx, 'pValue'] == 0) {
            dev[idx, 'pValue'] = 1e-4
        }
        df[i, j] = (sign(dev[idx, 'normalizedEnrichmentScore']) *
                    -log10(dev[idx, 'pValue']))
        pvals[i, j] = dev[idx, 'pValue']
    }
    mylim = max(abs(df))
    colnames(df)[1] = 'pan-developmental'

    print(sprintf('%.2f, %.2f', min(df), max(df)))

    # i checked and it comes out in the correct order
    pvals2 = matrix(p.adjust(pvals, method='fdr'), ncol=6, nrow=2)
    colnames(pvals2) = colnames(df)
    rownames(pvals2) = rownames(df)

    # figuring out rows and columns to highlight
    markme = which(round(pvals2, 2) <= .05, arr.ind = T) 

    quartz(width=5, height=3)
    corrplot(df, is.corr=F, col.lim=c(-mylim, mylim), tl.col='black',
            cl.length=3, cl.align.text='l', cl.offset=.2,
            col=colorRampPalette(c("#dc3220", "white","#005AB5"))(200)) -> p
    # make thicker borders for each significant square
    for (row in 1:nrow(markme)) {
        r = c(colnames(df)[markme[row, 'col']], rownames(df)[markme[row, 'row']],
            colnames(df)[markme[row, 'col']], rownames(df)[markme[row, 'row']])
        corrRect(p, namesMat = r) -> p
    }
    text(x=7.15, y=1.5, label=bquote(~log[10]~italic((p))), srt=270)

    p = recordPlot()
    return(p)
}

library(cowplot)
p1 = make_corrplot('main')
p2 = make_corrplot('WNH') 
sfig = plot_grid(p1, p2, rel_heights = c(1, 1), ncol=1,
                labels=c('A. Entire cohort',
                         'B. White non-Hispanics only'),
                label_x = 0, hjust=0, label_y=.9)
ggsave("~/data/post_mortem_adhd/figures/sfigure9.pdf", sfig, width=7.5, height=9,
       units="in")
       
```

## SFigure 10

```r
g1 = 'main'
g1_title = 'Entire cohort'
g2 = 'WNH'
g2_title = 'White non-Hispanics only'

tsize = 12
ysize = 10
xsize = 9
msize = .8

mydir = '~/data/post_mortem_adhd/results/'
r = 'ACC'
df = read.csv(sprintf('%s/GSEA_%s_%s_geneontology_Molecular_Function_noRedundant_10K.csv',
                      mydir, g1, r))
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore',
                         'FDR', 'pValue')]
df = df[round(df$FDR, 2) <= .05, ]
my_levels = df$description[sort(df$FDR, index.return=T, decreasing=T)$ix]
df$description = factor(df$description, levels=my_levels)
color_me = c('neurotransmitter receptor activity',
             'serotonin receptor activity', 'GABA receptor activity')
label_colors <- ifelse(levels(df$description) %in% color_me,"#DC3220",  "grey20")

df$Groups = g1_title
all_df = df

df = read.csv(sprintf('%s/GSEA_%s_%s_geneontology_Molecular_Function_noRedundant_10K.csv',
                      mydir, g2, r))
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore',
                            'FDR', 'pValue')]
df = df[df$description %in% my_levels, ]
df$description = factor(df$description, levels=my_levels)
df$Groups = g2_title
all_df = rbind(all_df, df)
all_df$Groups = factor(all_df$Groups)
all_df[all_df$pValue == 0, 'pValue'] = 1e-4

p <- ggplot(all_df,
            aes(y=-log10(pValue), x=description)) +
  geom_pointrange(aes(ymin=-log10(pValue), ymax=-log10(pValue), shape=Groups, colour=Groups),
                  size=msize) + coord_flip() +
  geom_hline(yintercept=-log10(.01), linetype="dotted",
                                color = "black", size=1) + 
  geom_hline(yintercept=-log10(.05), linetype="dashed",
                                color = "black", size=1) + 
    theme(axis.text.y = element_text(size=ysize, colour = label_colors),
          axis.title.y = element_blank(),
          plot.title = element_text(size = tsize),
          axis.text.x = element_text(size=xsize),
          axis.title.x = element_text(size=ysize)) +
    scale_colour_manual(values=c('#40B0A6', '#E1BE6A'))
p1 = p + ggtitle(r) + ylab(bquote(~-Log[10]~italic(P))) + ylim(0.25, 4.5)

r = 'Caudate'
df = read.csv(sprintf('%s/GSEA_%s_%s_geneontology_Molecular_Function_noRedundant_10K.csv',
                      mydir, g1, r))
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore',
                         'FDR', 'pValue')]
df = df[round(df$FDR, 2) <= .05, ]
my_levels = df$description[sort(df$FDR, index.return=T, decreasing=T)$ix]
df$description = factor(df$description, levels=my_levels)
color_me = c('fatty acid derivative binding', 'dynein light chain binding',
             'tau-protein kinase activity')
label_colors <- ifelse(levels(df$description) %in% color_me, "grey20", "#DC3220")
df$Groups = g1_title
all_df = df

df = read.csv(sprintf('%s/GSEA_%s_%s_geneontology_Molecular_Function_noRedundant_10K.csv',
                      mydir, g2, r))
df = df[order(df$FDR), c('description', 'normalizedEnrichmentScore',
                            'FDR', 'pValue')]
df = df[df$description %in% my_levels, ]
df$description = factor(df$description, levels=my_levels)
df$Groups = g2_title
all_df = rbind(all_df, df)
all_df$Groups = factor(all_df$Groups)
all_df[all_df$pValue == 0, 'pValue'] = 1e-4

p <- ggplot(all_df,
            aes(y=-log10(pValue), x=description)) +
  geom_pointrange(aes(ymin=-log10(pValue), ymax=-log10(pValue), shape=Groups, colour=Groups),
                  size=msize) + coord_flip() +
  geom_hline(yintercept=-log10(.01), linetype="dotted",
                                color = "black", size=1) + 
  geom_hline(yintercept=-log10(.05), linetype="dashed",
                                color = "black", size=1) + 
    theme(axis.text.y = element_text(size=ysize, colour = label_colors),
          axis.title.y = element_blank(),
          plot.title = element_text(size = tsize),
          axis.text.x = element_text(size=xsize),
          axis.title.x = element_text(size=ysize)) +
    scale_colour_manual(values=c('#40B0A6', '#E1BE6A'))
p2 = p + ggtitle(r) + ylab(bquote(~-Log[10]~italic(P))) + ylim(0.25, 4.5)

p = ggarrange(p2, p1, common.legend=T, ncol=2, nrow=1, legend='bottom')

ggsave("~/data/post_mortem_Adhd/figures/sfigure10.pdf", p, width=7.5, height=4,
       units="in")
```

## SFigure 11

```r
library(corrplot)

make_magma = function(g, cl_min=0, cl_max=8.5) {
    res = readRDS(sprintf('~/data/post_mortem_adhd/results/MAGMA_%s_res.rds', g))
    res = res[res$POP == g, c('REGION', 'DISORDER', 'P')]

    plot_mat = matrix(nrow=2, ncol=length(unique(res$DISORDER)))
    colnames(plot_mat) = unique(res$DISORDER)
    rownames(plot_mat) = unique(res$REGION)
    pvals = plot_mat
    for (i in 1:nrow(res)) {
        plot_mat[res[i, 'REGION'], res[i, 'DISORDER']] = -log10(res[i, 'P'])
        pvals[res[i, 'REGION'], res[i, 'DISORDER']] = res[i, 'P']
    }
    print(sprintf('%.2f, %.2f', min(plot_mat), max(plot_mat)))
    quartz()
    corrplot(plot_mat, is.corr=F, tl.col='black', p.mat = pvals,
            sig.level = .05, col.lim=c(cl_min, cl_max), cl.length=2,
            cl.align.text='l', cl.offset=.2, tl.cex = .8,
            col=colorRampPalette(c("white","#005AB5"))(200),
            mar = c(0, 0, 2, 0))
    text(x=10.25, y=1.5, label=bquote(~-log[10]~italic((p))), srt=270)
    magma = recordPlot()
    return(magma)
}

# adding some fake margin
library(cowplot)
p1 = make_magma('main')
p2 = make_magma('WNH')
sfig = plot_grid(p1, p2, rel_heights = c(1, 1), ncol=1,
                labels=c('A. Entire cohort',
                         'B. White non-Hispanics only'),
                label_x = 0, hjust=0, label_y=.95)
ggsave("~/data/post_mortem_adhd/figures/sfigure11.pdf", sfig, width=7.5, height=5,
       units="in")
```

## SFigure 12

```r
library(ggplot2)
library(ggpubr)

g1 = 'main'
g1_title = 'Entire Cohort'
g2 = 'WNH'
g2_title = 'White non-Hispanic'

# just to share axis
ymax = .5
ymin = -.15
leg_size = 9
cap_size = 10
tick_size = 9
title_size = 11

r = 'ACC'

corrs = readRDS(sprintf('~/data/post_mortem_adhd/results/disorders_corrs_%s.rds',
                        g1))
sources = unique(corrs$source)
dis_order = c('ASD', 'SCZ', 'BD', 'MDD', 'AAD', 'OCD')
col_labels = c('Autism Spectrum Disorder', 'Schizophrenia', 'Bipolar Disorder',
               'Major Depression Disorder', 'Alcohol Abuse or Dependence',
               'Obsessive Compulsive Disorder')
corrs$Disorders = factor(corrs$disorder, levels=dis_order)
# Bang Wong's color pallete
my_colors = c('#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2')
mycorrs = corrs[corrs$region == r, ]
mycorrs$id = sapply(1:nrow(mycorrs),
                  function(i) sprintf('%s [%d]',
                                      mycorrs[i, 'disorder'],
                                      which(sources == mycorrs[i, 'source'])))
df = data.frame()
for (st in unique(mycorrs$id)) {
    idx = mycorrs$id == st
    res = list(st=st, dis=as.character(unique(mycorrs[idx, 'disorder'])),
               val=median(mycorrs[idx, 'corr']),
               err=sd(mycorrs[idx, 'corr']), Groups=g1_title)
    df = rbind(df, res)
}
df$Disorders = factor(df$dis, level=dis_order)
all_df = df

corrs = readRDS(sprintf('~/data/post_mortem_adhd/results/disorders_corrs_%s.rds',
                        g2))
sources = unique(corrs$source)
mycorrs = corrs[corrs$region == r, ]
mycorrs$id = sapply(1:nrow(mycorrs),
                  function(i) sprintf('%s [%d]',
                                      mycorrs[i, 'disorder'],
                                      which(sources == mycorrs[i, 'source'])))

df = data.frame()
for (st in unique(mycorrs$id)) {
    idx = mycorrs$id == st
    res = list(st=st, dis=as.character(unique(mycorrs[idx, 'disorder'])),
            val=median(mycorrs[idx, 'corr']),
            err=sd(mycorrs[idx, 'corr']), Groups=g2_title)
    df = rbind(df, res)
}
df$Disorders = factor(df$dis, level=dis_order)
all_df = rbind(all_df, df)
all_df$Groups = factor(all_df$Groups)

mylevels = c()
ranks = all_df[all_df$Groups==g1_title, ]
for (d in dis_order) {
    # assumes IDs start with the disorder
    myranks = ranks[grepl(ranks$st, pattern=paste0('^',d)),]
    mylevels = c(mylevels, myranks$st[order(myranks$val, decreasing=T)])
}
# regroup levels based on the order we established before
all_df$xorder = factor(all_df$st, levels=mylevels)

p = ggplot(all_df, aes(y=val, x=xorder, color=Disorders)) + 
        geom_pointrange(aes(ymin=val-2*err, ymax=val+2*err, shape=Groups),
                        size=1) +
        ylim(ymin, ymax) + 
        geom_hline(yintercept=0, linetype="dotted",
                    color = "#dc3220", size=1) +
        # fake continuous axis to add vertical lines later
        geom_line(aes(x = as.numeric(xorder), y=0), size = 1, color="#dc3220",
                      alpha=0) + 
        # vertical lines separating disorders
        geom_vline(xintercept=c(4.5, 7.5, 10.5, 12.5),
                    linetype="dashed", color = "grey", size=1)
p1 = p + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=tick_size),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size=leg_size),
          legend.title = element_text(size=title_size),
          axis.text.y = element_text(size=tick_size),
          plot.title = element_text(size=title_size)) +
    scale_color_manual(breaks = dis_order,
                     values = my_colors,
                     labels = col_labels,
                     drop=FALSE) +
    ggtitle(r)


# repeat the same stuff with some hard coded modifications for Caudate, due to 
# different disorder datasets
r = 'Caudate'

corrs = readRDS(sprintf('~/data/post_mortem_adhd/results/disorders_corrs_%s.rds',
                        g1))
corrs$Disorders = factor(corrs$disorder, levels=dis_order)
mycorrs = corrs[corrs$region == r, ]
mycorrs$id = sapply(1:nrow(mycorrs),
                  function(i) sprintf('%s [%d]',
                                      mycorrs[i, 'disorder'],
                                      which(sources == mycorrs[i, 'source'])))
df = data.frame()
for (st in unique(mycorrs$id)) {
    idx = mycorrs$id == st
    res = list(st=st, dis=as.character(unique(mycorrs[idx, 'disorder'])),
               val=median(mycorrs[idx, 'corr']),
               err=sd(mycorrs[idx, 'corr']), Groups=g1_title)
    df = rbind(df, res)
}
df$Disorders = factor(df$dis, level=dis_order)
all_df = df

corrs = readRDS(sprintf('~/data/post_mortem_adhd/results/disorders_corrs_%s.rds',
                        g2))
mycorrs = corrs[corrs$region == r, ]
mycorrs$id = sapply(1:nrow(mycorrs),
                  function(i) sprintf('%s [%d]',
                                      mycorrs[i, 'disorder'],
                                      which(sources == mycorrs[i, 'source'])))

df = data.frame()
for (st in unique(mycorrs$id)) {
    idx = mycorrs$id == st
    res = list(st=st, dis=as.character(unique(mycorrs[idx, 'disorder'])),
            val=median(mycorrs[idx, 'corr']),
            err=sd(mycorrs[idx, 'corr']), Groups=g2_title)
    df = rbind(df, res)
}
df$Disorders = factor(df$dis, level=dis_order)
all_df = rbind(all_df, df)
all_df$Groups = factor(all_df$Groups)

mylevels = c()
ranks = all_df[all_df$Groups==g1_title, ]
for (d in dis_order) {
    # assumes IDs start with the disorder
    myranks = ranks[grepl(ranks$st, pattern=paste0('^',d)),]
    mylevels = c(mylevels, myranks$st[order(myranks$val, decreasing=T)])
}
# regroup levels based on the order we established before
all_df$xorder = factor(all_df$st, levels=mylevels)

p = ggplot(all_df, aes(y=val, x=xorder, color=Disorders)) + 
        geom_pointrange(aes(ymin=val-2*err, ymax=val+2*err, shape=Groups),
                        size=1) + ylim(ymin, ymax) + 
        geom_hline(yintercept=0, linetype="dotted",
                    color = "#dc3220", size=1) +
        # fake continuous axis to add vertical lines later
        geom_line(aes(x = as.numeric(xorder), y=0), size = 1, color="#dc3220",
                      alpha=0) + 
        # vertical lines separating disorders
        geom_vline(xintercept=c(2.5, 5.5, 8.5, 9.5, 10.5),
                    linetype="dashed", color = "grey", size=1)
p2 = p + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=tick_size),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size=leg_size),
          legend.title = element_text(size=title_size),
          axis.text.y = element_text(size=tick_size),
          plot.title = element_text(size=title_size)) +
    scale_color_manual(breaks = dis_order,
                     values = my_colors,
                     labels = col_labels,
                     drop=FALSE) +
    ggtitle(r)

p3 = ggarrange(p1, p2, common.legend = T, legend='right', nrow=2, ncol=1,
          legend.grob=get_legend(p2)) 
p4 = annotate_figure(p3, left = text_grob("Transcriptome correlation (rho)",
                                         rot = 90, size=cap_size))
ggsave("~/data/post_mortem_adhd/figures/sfigure12.pdf", p4, width=7, height=7,
       units="in")
```

## SFigure 13

It's a repeat of SFig9, but using the lessCov results.

```r
library(cowplot)
p1 = make_corrplot('main')
p2 = make_corrplot('lessCov') 
sfig = plot_grid(p1, p2, rel_heights = c(1, 1), ncol=1,
                labels=c('A. Entire cohort',
                         'B. Fewer covariates'),
                label_x = 0, hjust=0, label_y=.9)
ggsave("~/data/post_mortem_adhd/figures/sfigure13.pdf", sfig, width=7.5, height=9,
       units="in")
```

## SFigure 14

It's a repeat of SFig10, but using the lessCov results.

```r
g1 = 'main'
g1_title = 'Entire cohort'
g2 = 'lessCov'
g2_title = 'Fewer covariates'

# ...

ggsave("~/data/post_mortem_Adhd/figures/sfigure14.pdf", p, width=7.5, height=4,
       units="in")
```

## SFigure 15

It's a repeat of SFig11, but using the lessCov results.

```r
library(cowplot)
p1 = make_magma('main')
p2 = make_magma('lessCov')
sfig = plot_grid(p1, p2, rel_heights = c(1, 1), ncol=1,
                labels=c('A. Entire cohort',
                         'B. Fewer covariates'),
                label_x = 0, hjust=0, label_y=.95)
ggsave("~/data/post_mortem_adhd/figures/sfigure15.pdf", sfig, width=7.5, height=5,
       units="in")
```

## SFigure 16

It's a repeat of SFig12, but using the lessCov results.

```r
g1 = 'main'
g1_title = 'Entire Cohort'
g2 = 'lessCov'
g2_title = 'Fewer covariates'

# ...

ggsave("~/data/post_mortem_adhd/figures/sfigure16.pdf", p4, width=7, height=7,
       units="in")
```

## SFigure 17

It's a repeat of SFig9, but using the noMDD results.

```r
library(cowplot)
p1 = make_corrplot('main')
p2 = make_corrplot('noMDD') 
sfig = plot_grid(p1, p2, rel_heights = c(1, 1), ncol=1,
                labels=c('A. Entire cohort',
                         'B. No subjects with MDD'),
                label_x = 0, hjust=0, label_y=.9)
ggsave("~/data/post_mortem_adhd/figures/sfigure17.pdf", sfig, width=7.5, height=9,
       units="in")
```

## SFigure 18

It's a repeat of SFig10, but using the noMDD results.

```r
g1 = 'main'
g1_title = 'Entire cohort'
g2 = 'noMDD'
g2_title = 'No subjects with MDD'

# ...

ggsave("~/data/post_mortem_Adhd/figures/sfigure18.pdf", p, width=7.5, height=4,
       units="in")
```

## SFigure 19

It's a repeat of SFig11, but using the noMDD results.

```r
library(cowplot)
p1 = make_magma('main', cl_min=0, cl_max=8.6)
p2 = make_magma('noMDD', cl_min=0, cl_max=8.6)
sfig = plot_grid(p1, p2, rel_heights = c(1, 1), ncol=1,
                labels=c('A. Entire cohort',
                         'B. No subjects with MDD'),
                label_x = 0, hjust=0, label_y=.95)
ggsave("~/data/post_mortem_adhd/figures/sfigure19.pdf", sfig, width=7.5, height=5,
       units="in")
```

## SFigure 20

It's a repeat of SFig12, but using the noMDD results.

```r
g1 = 'main'
g1_title = 'Entire Cohort'
g2 = 'noMDD'
g2_title = 'No subjects with MDD'

# ...

ggsave("~/data/post_mortem_adhd/figures/sfigure20.pdf", p4, width=7, height=7,
       units="in")
```
