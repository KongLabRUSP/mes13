
> # Source: 
> # https://bioconductor.org/packages/release/bioc/html/DEGseq.html
> 
> # if (!requireNamespace("BiocManager",
> #                       .... [TRUNCATED] 

> require(ggplot2)

> require(DESeq2)

> require(DEGseq)

> require(knitr)

> require(ggdendro)

> # MES13 data----
> # Treatment legend----
> trt.names <- c("LG",
+                "HG",
+                "MIC_1.5uM",
+                "TIIA_5uM",
+ .... [TRUNCATED] 

> # Load data----
> # Question to Renyi: how was the data processed and annotated?
> # dt1 <- fread("data/Renyi_RNAseq_12292017/mes13_fpkm_Dec2017_dav ..." ... [TRUNCATED] 

> # dt1 <- fread("data/rna_seq/Renyi_12292017/mes13_featurecounts_Dec2017_david.csv")
> dt1
        Geneid                                                                                                 Chr
    1:    Xkr4                                                                                      chr1;chr1;chr1
    2:     Rp1                                                             chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1
    3:   Sox17 chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1
    4:  Mrpl15                               chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1
    5:  Lypla1                                                        chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1
   ---                                                                                                            
24417: Gm20816                                                                                 chrY;chrY;chrY;chrY
24418: Gm20867                                                                                 chrY;chrY;chrY;chrY
24419: Gm20806                                                                                 chrY;chrY;chrY;chrY
24420: Gm20854                                                                                 chrY;chrY;chrY;chrY
24421:   Erdr1                                                                                      chrY;chrY;chrY
                                                                                                                                                                 Start
    1:                                                                                                                                         3214482;3421702;3670552
    2:                                                                                                 4290846;4343507;4351910;4351910;4352202;4352202;4360200;4409170
    3: 4490928;4490928;4490928;4490928;4490928;4493100;4493100;4493772;4493772;4493772;4493772;4493772;4495136;4495136;4495136;4496291;4496291;4496291;4496291;4496291
    4:                                                 4773200;4773200;4773200;4777525;4777525;4777525;4782568;4782568;4782568;4783951;4783951;4785573;4785573;4785573
    5:                                                                                         4807893;4808455;4828584;4830268;4832311;4837001;4839387;4840956;4844963
   ---                                                                                                                                                                
24417:                                                                                                                             64301783;64561152;64784408;64786176
24418:                                                                                                                             64301785;64331715;64784408;64786176
24419:                                                                                                                             67339049;78835721;79958259;80756212
24420:                                                                                                                             83789754;83852592;83927454;85528525
24421:                                                                                                                                      90785442;90793296;90816349
                                                                                                                                                                   End
    1:                                                                                                                                         3216968;3421901;3671498
    2:                                                                                                 4293012;4350091;4352081;4352081;4352837;4352837;4360314;4409241
    3: 4492668;4492668;4492668;4492668;4492668;4493466;4493490;4493863;4493863;4493863;4493863;4493863;4495198;4495942;4495942;4497354;4497354;4497354;4497354;4497354
    4:                                                 4774516;4774516;4776801;4777648;4777648;4777648;4782733;4782733;4782733;4784105;4784105;4785726;4785726;4785726
    5:                                                                                         4807982;4808486;4828649;4830315;4832381;4837074;4839488;4841132;4846735
   ---                                                                                                                                                                
24417:                                                                                                                             64302222;64561152;64785064;64786299
24418:                                                                                                                             64302217;64331720;64785064;64786299
24419:                                                                                                                             67340047;78836719;79959257;80757210
24420:                                                                                                                             83790098;83852777;83927911;85529519
24421:                                                                                                                                      90785979;90793417;90816465
                                        Strand Length WJ1.dedup.bam WJ2.dedup.bam WJ3.dedup.bam WJ4.dedup.bam WJ5.dedup.bam WJ6.dedup.bam
    1:                                   -;-;-   3634             0             0             0             0             0             0
    2:                         -;-;-;-;-;-;-;-   9747             0             0             0             0             0             0
    3: -;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-   4095             0             0             0             0             1             0
    4:             -;-;-;-;-;-;-;-;-;-;-;-;-;-   4201           629           648           599           580           573           657
    5:                       +;+;+;+;+;+;+;+;+   2433          1143          1152          1032          1053          1170          1119
   ---                                                                                                                                   
24417:                                 -;-;-;-   1222             0             0             0             0             0             0
24418:                                 -;-;-;-   1220             0             0             0             0             0             0
24419:                                 -;-;-;-   3996             0             0             0             0             0             0
24420:                                 -;-;-;-   1984             0             0             0             0             0             0
24421:                                   +;+;+    777           428           390           395           456           362           459
       WJ7.dedup.bam
    1:             0
    2:             0
    3:             0
    4:           630
    5:          1167
   ---              
24417:             0
24418:             0
24419:             0
24420:             0
24421:           406

> # CHECK
> dt1[Geneid == "Tnfrsf25",]
     Geneid                                                                                  Chr
1: Tnfrsf25 chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4;chr4
                                                                                                                                                                       Start
1: 152115934;152116343;152116585;152116931;152116931;152117378;152117378;152117640;152118247;152118401;152118401;152118690;152118690;152119160;152119160;152119511;152119511
                                                                                                                                                                         End
1: 152116114;152116723;152116723;152117065;152117065;152117542;152117542;152117694;152118302;152118499;152118499;152118727;152118727;152119343;152119343;152120119;152120119
                              Strand Length WJ1.dedup.bam WJ2.dedup.bam WJ3.dedup.bam WJ4.dedup.bam WJ5.dedup.bam WJ6.dedup.bam
1: +;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+   1903             7             1             0             7             2             3
   WJ7.dedup.bam
1:             3

> dt1[Geneid == "Fcer1g",]
   Geneid                      Chr                                             Start                                               End
1: Fcer1g chr1;chr1;chr1;chr1;chr1 171229572;171230067;171230996;171231177;171234197 171229927;171230087;171231031;171231268;171234349
      Strand Length WJ1.dedup.bam WJ2.dedup.bam WJ3.dedup.bam WJ4.dedup.bam WJ5.dedup.bam WJ6.dedup.bam WJ7.dedup.bam
1: -;-;-;-;-    658             6             1             3             6             4             4             3

> # Keep only gene IDs and counts
> dt1 <- dt1[, c("Geneid",
+                "WJ1.dedup.bam",
+                "WJ2.dedup.bam",
+                "WJ3 ..." ... [TRUNCATED] 

> dt1
        Geneid WJ1.dedup.bam WJ2.dedup.bam WJ3.dedup.bam WJ4.dedup.bam WJ5.dedup.bam WJ6.dedup.bam WJ7.dedup.bam
    1:    Xkr4             0             0             0             0             0             0             0
    2:     Rp1             0             0             0             0             0             0             0
    3:   Sox17             0             0             0             0             1             0             0
    4:  Mrpl15           629           648           599           580           573           657           630
    5:  Lypla1          1143          1152          1032          1053          1170          1119          1167
   ---                                                                                                          
24417: Gm20816             0             0             0             0             0             0             0
24418: Gm20867             0             0             0             0             0             0             0
24419: Gm20806             0             0             0             0             0             0             0
24420: Gm20854             0             0             0             0             0             0             0
24421:   Erdr1           428           390           395           456           362           459           406

> colnames(dt1) <- c("GeneNames",
+                    trt.names)

> dt1
       GeneNames   LG   HG MIC_1.5uM TIIA_5uM FX_1uM Gen_10uM Ber_6uM
    1:      Xkr4    0    0         0        0      0        0       0
    2:       Rp1    0    0         0        0      0        0       0
    3:     Sox17    0    0         0        0      1        0       0
    4:    Mrpl15  629  648       599      580    573      657     630
    5:    Lypla1 1143 1152      1032     1053   1170     1119    1167
   ---                                                               
24417:   Gm20816    0    0         0        0      0        0       0
24418:   Gm20867    0    0         0        0      0        0       0
24419:   Gm20806    0    0         0        0      0        0       0
24420:   Gm20854    0    0         0        0      0        0       0
24421:     Erdr1  428  390       395      456    362      459     406

> # Remove genes with low counts----
> summary(dt1[, -1])
       LG                HG            MIC_1.5uM          TIIA_5uM           FX_1uM           Gen_10uM          Ber_6uM       
 Min.   :    0.0   Min.   :    0.0   Min.   :    0.0   Min.   :    0.0   Min.   :    0.0   Min.   :    0.0   Min.   :    0.0  
 1st Qu.:    0.0   1st Qu.:    0.0   1st Qu.:    0.0   1st Qu.:    0.0   1st Qu.:    0.0   1st Qu.:    0.0   1st Qu.:    0.0  
 Median :   14.0   Median :   15.0   Median :   15.0   Median :   14.0   Median :   16.0   Median :   17.0   Median :   16.0  
 Mean   :  386.6   Mean   :  397.8   Mean   :  386.3   Mean   :  398.1   Mean   :  394.5   Mean   :  402.1   Mean   :  405.9  
 3rd Qu.:  564.0   3rd Qu.:  575.0   3rd Qu.:  563.0   3rd Qu.:  585.0   3rd Qu.:  575.0   3rd Qu.:  581.0   3rd Qu.:  584.0  
 Max.   :10364.0   Max.   :10867.0   Max.   :10640.0   Max.   :11878.0   Max.   :11031.0   Max.   :11598.0   Max.   :11387.0  

> tmp <- rowSums(dt1[, -1])

> # Remove if total across 3 samples is no more than 10
> dt1 <- droplevels(subset(dt1,
+                          tmp > 20))

> dt1
       GeneNames   LG   HG MIC_1.5uM TIIA_5uM FX_1uM Gen_10uM Ber_6uM
    1:    Mrpl15  629  648       599      580    573      657     630
    2:    Lypla1 1143 1152      1032     1053   1170     1119    1167
    3:     Tcea1 1195 1226      1172     1120   1239     1221    1180
    4:     Rgs20   22   28        10       12     23       30      16
    5:   Atp6v1h 1056 1116      1111      976   1033     1091    1063
   ---                                                               
13950:     Kdm5d  284  292       314      390    287      322     296
13951:   Eif2s3y  356  381       342      421    346      337     332
13952:       Uty  284  325       257      390    294      301     304
13953:     Ddx3y  640  767       610      807    673      614     665
13954:     Erdr1  428  390       395      456    362      459     406

> # 13,954 genes left, down from 24,421 genes
> 
> # Leave the 3 treatments only----
> dt1 <- subset(dt1,
+               select = c("GeneNames",
+    .... [TRUNCATED] 

> dt1
       GeneNames   LG   HG TIIA_5uM
    1:    Mrpl15  629  648      580
    2:    Lypla1 1143 1152     1053
    3:     Tcea1 1195 1226     1120
    4:     Rgs20   22   28       12
    5:   Atp6v1h 1056 1116      976
   ---                             
13950:     Kdm5d  284  292      390
13951:   Eif2s3y  356  381      421
13952:       Uty  284  325      390
13953:     Ddx3y  640  767      807
13954:     Erdr1  428  390      456

> # Normalize data to fragments per million (FPM)----
> mat <- data.frame(sample = colnames(dt1)[-1],
+                   trt = colnames(dt1)[-1],
+   .... [TRUNCATED] 

> mat
    sample      trt repl
1       LG       LG    1
2       HG       HG    1
3 TIIA_5uM TIIA_5uM    1

> dtm <- as.matrix(dt1[, -1, with = FALSE])

> rownames(dtm) <- dt1$GeneNames

> head(dtm)
          LG   HG TIIA_5uM
Mrpl15   629  648      580
Lypla1  1143 1152     1053
Tcea1   1195 1226     1120
Rgs20     22   28       12
Atp6v1h 1056 1116      976
Rb1cc1  1796 1744     1516

> dds <- DESeqDataSetFromMatrix(countData = dtm, 
+                               colData = mat,
+                               ~ trt)

> dds <- estimateSizeFactors(dds)

> dds
class: DESeqDataSet 
dim: 13954 3 
metadata(1): version
assays(1): counts
rownames(13954): Mrpl15 Lypla1 ... Ddx3y Erdr1
rowData names(0):
colnames(3): LG HG TIIA_5uM
colData names(4): sample trt repl sizeFactor

> # Fragments per million (FPM) normalization----
> dt1.fpm <- data.table(GeneNames = dt1$GeneNames,
+                       fpm(dds,
+                .... [TRUNCATED] 

> colnames(dt1.fpm)[-1] <- paste(colnames(dt1.fpm)[-1],
+                                "fpm",
+                                sep = "_")

> dt1.fpm
       GeneNames     LG_fpm    HG_fpm TIIA_5uM_fpm
    1:    Mrpl15  66.636460  66.72155     59.67717
    2:    Lypla1 121.089783 118.61610    108.34494
    3:     Tcea1 126.598679 126.23553    115.23868
    4:     Rgs20   2.330687   2.88303      1.23470
    5:   Atp6v1h 111.872975 114.90934    100.42228
   ---                                            
13950:     Kdm5d  30.087050  30.06589     40.12775
13951:   Eif2s3y  37.714753  39.22980     43.31740
13952:       Uty  30.087050  33.46374     40.12775
13953:     Ddx3y  67.801803  78.97443     83.03358
13954:     Erdr1  45.342456  40.15649     46.91861

> # DEGseq----
> # a. (HG - LG)----
> DEGexp(geneExpMatrix1 = dt1,
+        geneCol1 = which(colnames(dt1) == "GeneNames"), 
+        expCol1 = which( .... [TRUNCATED] 
Please wait...
gene id column in geneExpMatrix1 for sample1:  1 
expression value column(s) in geneExpMatrix1: 3 
total number of reads uniquely mapped to genome obtained from sample1: 9712004 
gene id column in geneExpMatrix2 for sample2:  1 
expression value column(s) in geneExpMatrix2: 2 
total number of reads uniquely mapped to genome obtained from sample2: 9439277 

method to identify differentially expressed genes:  MARS 
qValue threshold (Storey et al. 2003): 0.5 
fold change: 1.231144 
log2 fold change: 0.3 
output directory: tmp/hg_lg 

Please wait ...
Identifying differentially expressed genes ...
Please wait patiently ...
output ...

Done ...
The results can be observed in directory:  tmp/hg_lg 

> hg_lg <- fread("tmp/hg_lg/output_score.txt")

> hg_lg
       GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized    z-score      p-value q-value(Benjamini et al. 1995)
    1:      Npr3    341    642        -0.9128016                   -0.9538942 -10.113830 4.797129e-24                   6.669928e-20
    2:     Acta2    241    449        -0.8976823                   -0.9387750  -8.345787 7.074134e-17                   4.917938e-13
    3:    Malat1   2498   1897         0.3970538                    0.3559611   8.141406 3.907150e-16                   1.810834e-12
    4: Trp53inp1   1882   1400         0.4268398                    0.3857471   7.618315 2.570082e-14                   8.933604e-11
    5:      Evpl    365    194         0.9118398                    0.8707472   6.962324 3.347045e-12                   7.756219e-09
   ---                                                                                                                              
13950:    Adam32      0      4                NA                           NA         NA           NA                             NA
13951:     Mpzl2      1      0                NA                           NA         NA           NA                             NA
13952:    Cyp1a1      2      0                NA                           NA         NA           NA                             NA
13953:       Nyx      0      1                NA                           NA         NA           NA                             NA
13954:      Phex      0      4                NA                           NA         NA           NA                             NA
       q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)
    1:                6.546674e-20                                         TRUE
    2:                4.827059e-13                                         TRUE
    3:                1.777371e-12                                         TRUE
    4:                8.768519e-11                                         TRUE
    5:                7.612891e-09                                         TRUE
   ---                                                                         
13950:                          NA                                        FALSE
13951:                          NA                                        FALSE
13952:                          NA                                        FALSE
13953:                          NA                                        FALSE
13954:                          NA                                        FALSE

> # CHECK:
> hg_lg[hg_lg$GeneNames == "Nmu", ]
   GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized  z-score     p-value q-value(Benjamini et al. 1995)
1:       Nmu      9      0          4.169925                     4.128832 2.947562 0.003202908                      0.1504501
   q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)
1:                   0.1476699                                         TRUE

> hg_lg[hg_lg$`Signature(q-value(Storey et al. 2003) < 0.5)`,]
     GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized    z-score      p-value q-value(Benjamini et al. 1995)
  1:      Npr3    341    642        -0.9128016                   -0.9538942 -10.113830 4.797129e-24                   6.669928e-20
  2:     Acta2    241    449        -0.8976823                   -0.9387750  -8.345787 7.074134e-17                   4.917938e-13
  3:    Malat1   2498   1897         0.3970538                    0.3559611   8.141406 3.907150e-16                   1.810834e-12
  4: Trp53inp1   1882   1400         0.4268398                    0.3857471   7.618315 2.570082e-14                   8.933604e-11
  5:      Evpl    365    194         0.9118398                    0.8707472   6.962324 3.347045e-12                   7.756219e-09
 ---                                                                                                                              
466:     Ltc4s     12      4         1.5849625                    1.5438698   1.991836 4.638903e-02                   5.068707e-01
467:     Dbil5     12      4         1.5849625                    1.5438698   1.991836 4.638903e-02                   5.068707e-01
468:   Col11a2     12      4         1.5849625                    1.5438698   1.991836 4.638903e-02                   5.068707e-01
469:   Igf2bp3    156    188        -0.2691866                   -0.3102793  -1.990260 4.656233e-02                   5.069715e-01
470:      Ydjc     44     63        -0.5178483                   -0.5589410  -1.987916 4.682104e-02                   5.093894e-01
     q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)
  1:                6.546674e-20                                         TRUE
  2:                4.827059e-13                                         TRUE
  3:                1.777371e-12                                         TRUE
  4:                8.768519e-11                                         TRUE
  5:                7.612891e-09                                         TRUE
 ---                                                                         
466:                4.967876e-01                                         TRUE
467:                4.967876e-01                                         TRUE
468:                4.967876e-01                                         TRUE
469:                4.976031e-01                                         TRUE
470:                4.999764e-01                                         TRUE

> # Write as CSV----
> write.csv(hg_lg,
+           file = "tmp/hg_lg/MES13_RNA_DEGseq_HG-LG.csv",
+           row.names = FALSE)

> # MA Plot----
> hg_lg <- merge(hg_lg,
+                dt1.fpm,
+                by = "GeneNames")

> hg_lg[, mu := (log2(HG_fpm + 1) + log2(LG_fpm + 1))/2]
           GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized     z-score     p-value q-value(Benjamini et al. 1995)
    1: 0610007P14Rik    456    494      -0.115477217                  -0.15656987 -1.67195427 0.094533341                      0.6565392
    2: 0610009B22Rik    271    255       0.087795604                   0.04670295  0.37115911 0.710519028                      0.9707238
    3: 0610009L18Rik      4      6      -0.584962501                  -0.62605516 -0.67920051 0.497010826                      0.9328346
    4: 0610009O20Rik    910    909       0.001586251                  -0.03950640 -0.58401222 0.559212093                      0.9457834
    5: 0610010B08Rik      4      4       0.000000000                  -0.04109265 -0.04028346 0.967867140                      0.9987550
   ---                                                                                                                                  
13950:        Zyg11b   1403   1311       0.097847323                   0.05675467  1.02454238 0.305579166                      0.8574718
13951:           Zyx   1546   1664      -0.106115114                  -0.14720777 -2.88994979 0.003853034                      0.1633310
13952:         Zzef1   1382   1250       0.144829521                   0.10373687  1.84352882 0.065251838                      0.5782419
13953:          Zzz3   2038   2029       0.006385186                  -0.03470747 -0.76722563 0.442947377                      0.9151174
13954:         l7Rn6    345    343       0.008387785                  -0.03270487 -0.29732471 0.766218615                      0.9767584
       q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)      LG_fpm      HG_fpm TIIA_5uM_fpm        mu
    1:                   0.6392561                                        FALSE  52.3345167  46.9522047   44.0376382 5.6602614
    2:                   0.9527641                                        FALSE  27.0147809  27.9036129   27.3691864 4.8306471
    3:                   0.9124584                                        FALSE   0.6356419   0.4118614    0.7202417 0.6037277
    4:                   0.9281658                                        FALSE  96.2997484  93.6984787  100.9367361 6.5848168
    5:                   0.9767325                                        FALSE   0.4237613   0.4118614    0.4115667 0.5036529
   ---                                                                                                                        
13950:                   0.8401522                                        FALSE 138.8877559 144.4604018  130.7753227 7.1563043
13951:                   0.1603128                                        FALSE 176.2846879 159.1844484  180.1633282 7.3967572
13952:                   0.5673136                                        FALSE 132.4253966 142.2981292  163.4948764 7.1113827
13953:                   0.8977791                                        FALSE 214.9529037 209.8434062  216.2783072 7.7373005
13954:                   0.9574830                                        FALSE  36.3375288  35.5230496   35.6005206 5.2066449

> hg_lg[, diff := log2(HG_fpm + 1) - log2(LG_fpm + 1)]
           GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized     z-score     p-value q-value(Benjamini et al. 1995)
    1: 0610007P14Rik    456    494      -0.115477217                  -0.15656987 -1.67195427 0.094533341                      0.6565392
    2: 0610009B22Rik    271    255       0.087795604                   0.04670295  0.37115911 0.710519028                      0.9707238
    3: 0610009L18Rik      4      6      -0.584962501                  -0.62605516 -0.67920051 0.497010826                      0.9328346
    4: 0610009O20Rik    910    909       0.001586251                  -0.03950640 -0.58401222 0.559212093                      0.9457834
    5: 0610010B08Rik      4      4       0.000000000                  -0.04109265 -0.04028346 0.967867140                      0.9987550
   ---                                                                                                                                  
13950:        Zyg11b   1403   1311       0.097847323                   0.05675467  1.02454238 0.305579166                      0.8574718
13951:           Zyx   1546   1664      -0.106115114                  -0.14720777 -2.88994979 0.003853034                      0.1633310
13952:         Zzef1   1382   1250       0.144829521                   0.10373687  1.84352882 0.065251838                      0.5782419
13953:          Zzz3   2038   2029       0.006385186                  -0.03470747 -0.76722563 0.442947377                      0.9151174
13954:         l7Rn6    345    343       0.008387785                  -0.03270487 -0.29732471 0.766218615                      0.9767584
       q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)      LG_fpm      HG_fpm TIIA_5uM_fpm        mu        diff
    1:                   0.6392561                                        FALSE  52.3345167  46.9522047   44.0376382 5.6602614 -0.15347236
    2:                   0.9527641                                        FALSE  27.0147809  27.9036129   27.3691864 4.8306471  0.04506163
    3:                   0.9124584                                        FALSE   0.6356419   0.4118614    0.7202417 0.6037277 -0.21225841
    4:                   0.9281658                                        FALSE  96.2997484  93.6984787  100.9367361 6.5848168 -0.03909482
    5:                   0.9767325                                        FALSE   0.4237613   0.4118614    0.4115667 0.5036529 -0.01210875
   ---                                                                                                                                    
13950:                   0.8401522                                        FALSE 138.8877559 144.4604018  130.7753227 7.1563043  0.05635677
13951:                   0.1603128                                        FALSE 176.2846879 159.1844484  180.1633282 7.3967572 -0.14633385
13952:                   0.5673136                                        FALSE 132.4253966 142.2981292  163.4948764 7.1113827  0.10298648
13953:                   0.8977791                                        FALSE 214.9529037 209.8434062  216.2783072 7.7373005 -0.03454481
13954:                   0.9574830                                        FALSE  36.3375288  35.5230496   35.6005206 5.2066449 -0.03181921

> tiff(filename = "tmp/hg_lg/MES13_RNA_DEGseq_HG-LG_maplot.tiff",
+      height = 6,
+      width = 6,
+      units = 'in',
+      res = 300,
+      c .... [TRUNCATED] 

> plot(hg_lg$diff ~ hg_lg$mu,
+      pch = ".",
+      xlab = "Mean",
+      ylab = "Difference",
+      main = "MES13 FPM-Normalized Gene Expressions ..." ... [TRUNCATED] 

> points(hg_lg$diff[hg_lg$`q-value(Storey et al. 2003)` < 0.5 & hg_lg$diff > 0] ~ hg_lg$mu[hg_lg$`q-value(Storey et al. 2003)` < 0.5 & hg_lg$diff > 0] .... [TRUNCATED] 

> points(hg_lg$diff[hg_lg$`q-value(Storey et al. 2003)` < 0.5 & hg_lg$diff < 0] ~ hg_lg$mu[hg_lg$`q-value(Storey et al. 2003)` < 0.5 & hg_lg$diff < 0] .... [TRUNCATED] 

> abline(h = c(-0.3, 0.3),
+        lty = 2)

> graphics.off()

> # b. (TIIA - HG)----
> DEGexp(geneExpMatrix1 = dt1,
+        geneCol1 = which(colnames(dt1) == "GeneNames"), 
+        expCol1 = which(colnames(dt1) .... [TRUNCATED] 
Please wait...
gene id column in geneExpMatrix1 for sample1:  1 
expression value column(s) in geneExpMatrix1: 4 
total number of reads uniquely mapped to genome obtained from sample1: 9718959 
gene id column in geneExpMatrix2 for sample2:  1 
expression value column(s) in geneExpMatrix2: 3 
total number of reads uniquely mapped to genome obtained from sample2: 9712004 

method to identify differentially expressed genes:  MARS 
qValue threshold (Storey et al. 2003): 0.5 
fold change: 1.231144 
log2 fold change: 0.3 
output directory: tmp/tiia_hg 

Please wait ...
Identifying differentially expressed genes ...
Please wait patiently ...
output ...

Done ...
The results can be observed in directory:  tmp/tiia_hg 

> tiia_hg <- fread("tmp/tiia_hg/output_score.txt")

> tiia_hg
       GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized   z-score      p-value q-value(Benjamini et al. 1995)
    1:        Cp   1160   2145        -0.8868528                   -0.8878856 -17.28519 6.082386e-67                   8.472156e-63
    2:      Scd2   2732   4059        -0.5711669                   -0.5721996 -16.18667 6.263044e-59                   4.361897e-55
    3:        C3     26    290        -3.4794694                   -3.4805022 -15.89634 6.717581e-57                   3.118973e-53
    4:      Lrp1   4388   3180         0.4645368                    0.4635040  13.88740 7.552944e-44                   2.630124e-40
    5:    Il1rl1    299    735        -1.2975988                   -1.2986315 -13.78150 3.293503e-43                   9.175042e-40
   ---                                                                                                                             
13950:    Trim66      0      3                NA                           NA        NA           NA                             NA
13951:  Tnfrsf26      0      4                NA                           NA        NA           NA                             NA
13952:   Cbfa2t3      0      4                NA                           NA        NA           NA                             NA
13953:   Gm14827      2      0                NA                           NA        NA           NA                             NA
13954:      Phex      3      0                NA                           NA        NA           NA                             NA
       q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)
    1:                6.039997e-63                                         TRUE
    2:                3.109698e-55                                         TRUE
    3:                2.223588e-53                                         TRUE
    4:                1.875077e-40                                         TRUE
    5:                6.541100e-40                                         TRUE
   ---                                                                         
13950:                          NA                                        FALSE
13951:                          NA                                        FALSE
13952:                          NA                                        FALSE
13953:                          NA                                        FALSE
13954:                          NA                                        FALSE

> #CHECK:
> tiia_hg[tiia_hg$GeneNames == "Nmu"]
   GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized   z-score   p-value q-value(Benjamini et al. 1995)
1:       Nmu      5      9        -0.8479969                   -0.8490297 -1.077796 0.2811247                      0.5862833
   q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)
1:                   0.4176623                                         TRUE

> tiia_hg[tiia_hg$`Signature(q-value(Storey et al. 2003) < 0.5)`, ]
          GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized     z-score      p-value q-value(Benjamini et al. 1995)
   1:            Cp   1160   2145        -0.8868528                   -0.8878856 -17.2851869 6.082386e-67                   8.472156e-63
   2:          Scd2   2732   4059        -0.5711669                   -0.5721996 -16.1866670 6.263044e-59                   4.361897e-55
   3:            C3     26    290        -3.4794694                   -3.4805022 -15.8963406 6.717581e-57                   3.118973e-53
   4:          Lrp1   4388   3180         0.4645368                    0.4635040  13.8873982 7.552944e-44                   2.630124e-40
   5:        Il1rl1    299    735        -1.2975988                   -1.2986315 -13.7815014 3.293503e-43                   9.175042e-40
  ---                                                                                                                                   
2509:         Naprt     14     10         0.4854268                    0.4843940   0.8166605 4.141225e-01                   6.973298e-01
2510: A230056P14Rik     14     10         0.4854268                    0.4843940   0.8166605 4.141225e-01                   6.973298e-01
2511:  LOC100504703     30     24         0.3219281                    0.3208953   0.8147185 4.152335e-01                   6.981034e-01
2512: 4931440P22Rik     30     24         0.3219281                    0.3208953   0.8147185 4.152335e-01                   6.981034e-01
2513:         Aldob     30     24         0.3219281                    0.3208953   0.8147185 4.152335e-01                   6.981034e-01
      q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)
   1:                6.039997e-63                                         TRUE
   2:                3.109698e-55                                         TRUE
   3:                2.223588e-53                                         TRUE
   4:                1.875077e-40                                         TRUE
   5:                6.541100e-40                                         TRUE
  ---                                                                         
2509:                4.970525e-01                                         TRUE
2510:                4.970525e-01                                         TRUE
2511:                4.976341e-01                                         TRUE
2512:                4.976341e-01                                         TRUE
2513:                4.976341e-01                                         TRUE

> # Write as CSV----
> write.csv(tiia_hg,
+           file = "tmp/tiia_hg/MES13_RNA_DEGseq_TIIA-HG.csv",
+           row.names = FALSE)

> # MA Plot----
> tiia_hg <- merge(tiia_hg,
+                  dt1.fpm,
+                  by = "GeneNames")

> tiia_hg[, mu := (log2(TIIA_5uM_fpm + 1) + log2(HG_fpm + 1))/2]
           GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized     z-score      p-value q-value(Benjamini et al. 1995)
    1: 0610007P14Rik    428    456      -0.091423028                 -0.092455807 -0.95248182 0.3408526801                    0.641412724
    2: 0610009B22Rik    266    271      -0.026866606                 -0.027899385 -0.22406490 0.8227068010                    0.928795837
    3: 0610009L18Rik      7      4       0.807354922                  0.806322143  0.90909379 0.3633006218                    0.658736574
    4: 0610009O20Rik    981    910       0.108386591                  0.107353812  1.61743437 0.1057845522                    0.339197290
    5: 0610010B08Rik      4      4       0.000000000                 -0.001032779 -0.00101239 0.9991922296                    1.000000000
   ---                                                                                                                                   
13950:        Zyg11b   1271   1403      -0.142550979                 -0.143583758 -2.57185480 0.0101155307                    0.079245910
13951:           Zyx   1751   1546       0.179638765                  0.178605986  3.55113491 0.0003835737                    0.007048547
13952:         Zzef1   1589   1382       0.201361509                  0.200328730  3.78002632 0.0001568118                    0.003450602
13953:          Zzz3   2102   2038       0.044608618                  0.043575839  0.97176605 0.3311669499                    0.631980332
13954:         l7Rn6    346    345       0.004175676                  0.003142897  0.02863335 0.9771570139                    1.000000000
       q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)      LG_fpm      HG_fpm TIIA_5uM_fpm        mu
    1:                 0.457278027                                        FALSE  52.3345167  46.9522047   44.0376382 5.5382923
    2:                 0.662160122                                        FALSE  27.0147809  27.9036129   27.3691864 4.8397154
    3:                 0.469327048                                         TRUE   0.6356419   0.4118614    0.7202417 0.6401049
    4:                 0.241731422                                        FALSE  96.2997484  93.6984787  100.9367361 6.6183998
    5:                 0.712670603                                        FALSE   0.4237613   0.4118614    0.4115667 0.4974479
   ---                                                                                                                        
13950:                 0.056496250                                        FALSE 138.8877559 144.4604018  130.7753227 7.1132095
13951:                 0.005024487                                        FALSE 176.2846879 159.1844484  180.1633282 7.4123687
13952:                 0.002458498                                        FALSE 132.4253966 142.2981292  163.4948764 7.2623874
13953:                 0.450553456                                        FALSE 214.9529037 209.8434062  216.2783072 7.7417142
13954:                 0.712670603                                        FALSE  36.3375288  35.5230496   35.6005206 5.1922638

> tiia_hg[, diff := log2(TIIA_5uM_fpm + 1) - log2(HG_fpm + 1)]
           GeneNames value1 value2 log2(Fold_change) log2(Fold_change) normalized     z-score      p-value q-value(Benjamini et al. 1995)
    1: 0610007P14Rik    428    456      -0.091423028                 -0.092455807 -0.95248182 0.3408526801                    0.641412724
    2: 0610009B22Rik    266    271      -0.026866606                 -0.027899385 -0.22406490 0.8227068010                    0.928795837
    3: 0610009L18Rik      7      4       0.807354922                  0.806322143  0.90909379 0.3633006218                    0.658736574
    4: 0610009O20Rik    981    910       0.108386591                  0.107353812  1.61743437 0.1057845522                    0.339197290
    5: 0610010B08Rik      4      4       0.000000000                 -0.001032779 -0.00101239 0.9991922296                    1.000000000
   ---                                                                                                                                   
13950:        Zyg11b   1271   1403      -0.142550979                 -0.143583758 -2.57185480 0.0101155307                    0.079245910
13951:           Zyx   1751   1546       0.179638765                  0.178605986  3.55113491 0.0003835737                    0.007048547
13952:         Zzef1   1589   1382       0.201361509                  0.200328730  3.78002632 0.0001568118                    0.003450602
13953:          Zzz3   2102   2038       0.044608618                  0.043575839  0.97176605 0.3311669499                    0.631980332
13954:         l7Rn6    346    345       0.004175676                  0.003142897  0.02863335 0.9771570139                    1.000000000
       q-value(Storey et al. 2003) Signature(q-value(Storey et al. 2003) < 0.5)      LG_fpm      HG_fpm TIIA_5uM_fpm        mu
    1:                 0.457278027                                        FALSE  52.3345167  46.9522047   44.0376382 5.5382923
    2:                 0.662160122                                        FALSE  27.0147809  27.9036129   27.3691864 4.8397154
    3:                 0.469327048                                         TRUE   0.6356419   0.4118614    0.7202417 0.6401049
    4:                 0.241731422                                        FALSE  96.2997484  93.6984787  100.9367361 6.6183998
    5:                 0.712670603                                        FALSE   0.4237613   0.4118614    0.4115667 0.4974479
   ---                                                                                                                        
13950:                 0.056496250                                        FALSE 138.8877559 144.4604018  130.7753227 7.1132095
13951:                 0.005024487                                        FALSE 176.2846879 159.1844484  180.1633282 7.4123687
13952:                 0.002458498                                        FALSE 132.4253966 142.2981292  163.4948764 7.2623874
13953:                 0.450553456                                        FALSE 214.9529037 209.8434062  216.2783072 7.7417142
13954:                 0.712670603                                        FALSE  36.3375288  35.5230496   35.6005206 5.1922638
                diff
    1: -0.0904659748
    2: -0.0269250642
    3:  0.2850128070
    4:  0.1062609113
    5: -0.0003012009
   ---              
13950: -0.1425462408
13951:  0.1775568588
13952:  0.1990228730
13953:  0.0433722406
13954:  0.0030569357

> tiff(filename = "tmp/tiia_hg/MES13_RNA_DEGseq_TIIA-HG_maplot.tiff",
+      height = 6,
+      width = 6,
+      units = 'in',
+      res = 300,
+    .... [TRUNCATED] 

> plot(tiia_hg$diff ~ tiia_hg$mu,
+      pch = ".",
+      xlab = "Mean",
+      ylab = "Difference",
+      main = "MES13 FPM-Normalized Gene Express ..." ... [TRUNCATED] 

> points(tiia_hg$diff[tiia_hg$`q-value(Storey et al. 2003)` < 0.5 & tiia_hg$diff > 0] ~ tiia_hg$mu[tiia_hg$`q-value(Storey et al. 2003)` < 0.5 & tiia_ .... [TRUNCATED] 

> points(tiia_hg$diff[tiia_hg$`q-value(Storey et al. 2003)` < 0.5 & tiia_hg$diff < 0] ~ tiia_hg$mu[tiia_hg$`q-value(Storey et al. 2003)` < 0.5 & tiia_ .... [TRUNCATED] 

> abline(h = c(-0.3, 0.3),
+        lty = 2)

> graphics.off()

> # Venn diagram----
> g1 <- hg_lg[`q-value(Storey et al. 2003)` < 0.5 & 
+               `log2(Fold_change) normalized` > 0.3,]$GeneNames

> # 263 genes
> g2 <- hg_lg[`q-value(Storey et al. 2003)` < 0.5 & 
+               `log2(Fold_change) normalized` < -0.3,]$GeneNames

> # 207 genes
> 
> g3 <- tiia_hg[`q-value(Storey et al. 2003)` < 0.5 & 
+                 `log2(Fold_change) normalized` > 0.3,]$GeneNames

> # 1,120 genes
> g4 <- tiia_hg[`q-value(Storey et al. 2003)` < 0.5 & 
+                 `log2(Fold_change) normalized` < -0.3,]$GeneNames

> # 1,393 genes
> 
> up.dn <- g1[g1 %in% g4]

> # 124 genes
> dn.up <- g2[g2 %in% g3]

> # 89 genes
> 
> # Heatmap----
> ll <- unique(c(up.dn,
+                dn.up))

> t1 <- merge(hg_lg[hg_lg$GeneNames %in% ll, 
+                   c("GeneNames",
+                     "log2(Fold_change) normalized")],
+             .... [TRUNCATED] 

> colnames(t1) <- c("Gene",
+                   "HG-LG",
+                   "TIIA-HG")

> t1 <- t1[order(t1$`HG-LG`,
+                decreasing = TRUE), ]

> t1
        Gene     HG-LG    TIIA-HG
  1:     Nmu  4.128832 -0.8490297
  2: Fam159b  3.958907 -0.5416012
  3:   Them5  3.958907 -4.0010328
  4:  Hsd3b2  3.659347 -2.1165100
  5:  Cd300a  3.543870 -3.5859953
 ---                             
209:  Kcnip2 -3.041093  1.9989672
210:   Tssk2 -3.041093  2.3208953
211:   Grin1 -3.500524  1.9989672
212:  Sh3bgr -3.626055  3.9989672
213:    Lyz2 -3.848448  3.3208953

> t1[t1$Gene == "Nmu", ]
   Gene    HG-LG    TIIA-HG
1:  Nmu 4.128832 -0.8490297

> # Gene    HG-LG      TIIA-HG
> # Nmu     4.128832   -0.8490297
> write.csv(t1,
+           file = "tmp/mes13_tiia_rnaseq_degseq_genes_q-0.5_log2-0.3 ..." ... [TRUNCATED] 

> ll <- melt.data.table(data = t1,
+                       id.vars = 1,
+                       measure.vars = 2:3,
+                       variable.n .... [TRUNCATED] 

> ll$Comparison <- factor(ll$Comparison,
+                         levels = c("TIIA-HG", 
+                                    "HG-LG"))

> lvls <- ll[ll$Comparison == "HG-LG", ]

> ll$Gene <- factor(ll$Gene,
+                   levels = lvls$Gene[order(lvls$`Gene Expression Diff`)])

> # Add dendrogram----
> # Sources: 
> # https://stackoverflow.com/questions/43794870/plotting-a-clustered-heatmap-with-dendrograms-using-rs-plotly
>  .... [TRUNCATED] 

> rownames(dt.dndr) <- dt.dndr$Gene

> dt.dndr <- dt.dndr[, -1]

> # Compute distances between genes----
> sampleDists <- dist(dt.dndr)

> as.matrix(sampleDists)
                    Nmu   Fam159b      Them5    Hsd3b2     Cd300a   Fam196a    Fbxl13       C4b    Entpd1     Scn1a   Snora5c    Zfp345
Nmu           0.0000000 0.3512646  3.1565801 1.3516371  2.7987786 0.6043890 0.9409035 1.2285618 1.7004623 1.7004623 0.9714471 1.7004623
Fam159b       0.3512646 0.0000000  3.4594316 1.6031451  3.0725546 0.6191394 1.1238395 1.3751901 1.9060493 1.9060493 1.0345552 1.9060493
Them5         3.1565801 3.4594316  0.0000000 1.9081830  0.5869517 3.0285733 2.4504412 2.3624365 1.8098914 1.8098914 2.7625804 1.8098914
Hsd3b2        1.3516371 1.6031451  1.9081830 0.0000000  1.4740156 1.1214385 0.5429372 0.5352963 0.4311531 0.4311531 0.8791992 0.4311531
                 Kcnj15     Obscn     Cd59a   Mkln1os    Nckap5 5330417C22Rik     Akap5     Sox21      Fgl2   Gm19589 4933433G19Rik
Nmu           1.0505416 2.5281120 1.1822009 1.3631761 1.4130708     2.1660700 1.5849625 1.6384023 1.9948516  3.996029     1.8542377
Fam159b       1.0416882 2.7572705 1.0094093 1.2221200 1.3737854     2.3052395 1.4480481 1.5888009 2.0378033  4.210800     1.7398284
Them5         2.9492964 1.1739033 3.4691795 3.4081351 3.0228308     2.0013605 3.4550622 3.0811218 2.6251157  1.607616     3.4374300
Hsd3b2        1.0842116 1.1803838 1.5989840 1.5839143 1.2630126     1.0366953 1.6884300 1.4039046 1.2641938  2.652211     1.7732974
                 Gm15708    Dmrta2     Phf11a    Slc5a1   Alox12b D630023F18Rik     Etnk2      Plb1     Ifi44    Gm3704   Slc35f3
Nmu           1.96844831 2.1314592 2.00000000 2.3973931 2.0840942     2.3268982 2.5499226 2.4402920 2.8023320 2.4839290 2.5258273
Fam159b       1.82042782 2.1071150 1.85571733 2.4503015 1.9257228     2.2004987 2.5322797 2.2661293 2.8239775 2.3703370 2.4474735
Them5         3.65057161 3.0301123 3.64476309 2.5881169 3.7576008     3.6920343 3.1008684 4.0353256 2.9068593 3.6889632 3.4871448
Hsd3b2        1.97904027 1.6198522 1.98720448 1.5314829 2.1033400     2.1623676 1.9179962 2.4447217 1.9861023 2.2360680 2.1314592
                Slc6a12     Aknad1      Catip   Col11a2     Dbil5     Ltc4s       Otog      Thbs2   Mab21l3      Ptges     Sprr1b
Nmu           2.5398510 2.62804907 2.58592068 2.8300428 2.5983106 2.5894277 2.58563331 2.63568145 2.6385246 2.65539142 2.66869900
Fam159b       2.3391733 2.53829256 2.42664248 2.8217631 2.4154455 2.4583498 2.44266209 2.55109641 2.5319926 2.50618918 2.50761130
Them5         4.2580813 3.60617183 4.02693668 3.1356668 4.1826890 3.8512863 3.92424922 3.57609322 3.7197632 4.00265848 4.08647403
Hsd3b2        2.6375384 2.25941675 2.50302597 2.1186266 2.6110762 2.3915546 2.43637548 2.24546499 2.3341398 2.51951846 2.57939570
                   Esam   Fam110c   Slc6a18     Ccr10    Rimbp3    Nipal1    Lrrc15     Itpka   Ptgs2os    Iigp1       Liph     Pkd1l3
Nmu           2.7591389 2.7218207 2.7147522 2.7549458 2.7674162 2.8901551 2.7955489 2.8193700 2.8611762 4.619445 2.96839536 2.96935083
Fam159b       2.6862308 2.5269476 2.5817756 2.6011238 2.5853437 2.8302338 2.6614430 2.6521852 2.6597343 4.764899 2.81259484 2.80764767
Them5         3.5561435 4.3375321 3.9312196 4.0902880 4.2830494 3.5290061 3.9839931 4.2174084 4.4733327 2.749155 4.23360430 4.27461318
Hsd3b2        2.3028777 2.7662273 2.5032792 2.6220836 2.7505524 2.3642123 2.5759735 2.7322958 2.9144536 3.362547 2.81322095 2.83864552
                   Gypc     Bambi     Efr3b Hist1h2ba       Gng2   Pacsin1     Inpp5d       Rgcc      Grem1    Pcdhgc4    Pcdhgb6
Nmu           2.9884006 3.2435740 3.0549081 3.0217474 3.12916912 3.1130120 3.14841198 3.18048160 3.22567104 3.24677371 3.23623552
Fam159b       2.8147867 3.2186714 2.9532371 2.8848161 3.01075690 2.9476564 2.96529455 3.05299013 3.11757791 3.13565510 3.05460094
Them5         4.3685627 3.4466840 3.9093677 4.1365549 4.07155638 4.3928920 4.54204314 4.16772596 4.05568781 4.09051071 4.59201673
Hsd3b2        2.9053496 2.5437253 2.6691290 2.7818732 2.80082567 2.9796073 3.08660177 2.88263283 2.84476925 2.87552750 3.15821643
                  Klra2       Evpl       Hey2 A630066F11Rik     Slc2a4      Tdrd5      Kcnj2      Esrp1     Thsd7a       Map7       Hic1
Nmu           3.2757124 3.27424769 3.29784589    3.31934028 3.30755191 3.32777174 3.32815374 3.33400253 3.33861095 3.34980850 3.35090399
Fam159b       3.0752501 3.08820977 3.13081939    3.13761520 3.14987157 3.17715828 3.16996353 3.17299538 3.17372853 3.17767485 3.18003597
Them5         4.7561185 4.65035372 4.52892617    4.65071378 4.46765802 4.42978431 4.48528725 4.50980293 4.54123661 4.60196788 4.59346852
Hsd3b2        3.2759896 3.21123912 3.15064914    3.23289001 3.11988238 3.10842381 3.14044142 3.15752830 3.17803059 3.21893527 3.21450903
                  Epcam    Gm15446  Rasgef1a     Alcam       Has3       Oas3    Fam212a      Grhl3      Adcy7        Fv1      Mdga1
Nmu           3.4720085 3.40534292 3.4439900 3.4943777 3.47276438 3.50107820 3.52246006 3.52697998 3.54209235 3.56944292 3.56219132
Fam159b       3.4049084 3.22679840 3.2453905 3.3824863 3.29420852 3.32235790 3.33544565 3.34835543 3.35230932 3.36648842 3.36952561
Them5         3.8978232 4.68808310 4.8648478 4.2524035 4.73618995 4.75772952 4.83547454 4.77566394 4.87071874 4.99090914 4.90734027
Hsd3b2        2.9063720 3.29643246 3.4193474 3.1005331 3.35736696 3.38372234 3.43927882 3.40680364 3.46930144 3.55293715 3.50033882
                  Clstn3     Dnah7b    Cyfip2       Ntm    Pip5k1b      Ptpn7     Celsr1        Sp6      Ism1       Shc4      Pax9
Nmu           3.55118822 3.57558168 3.5559995 3.5642673 3.58817948 3.61261558 3.61171319 3.63992643 3.6725212 3.61520114 3.6092634
Fam159b       3.36865488 3.38232677 3.3900345 3.4145366 3.40439895 3.41259511 3.42370983 3.43810065 3.5819252 3.43231585 3.4685739
Them5         4.82264494 4.92172516 4.7011866 4.5845967 4.85909825 5.00094662 4.90852658 5.03537958 4.2051593 4.87211661 4.5473442
Hsd3b2        3.44578929 3.51515435 3.3787815 3.3179113 3.48491064 3.57916581 3.52490408 3.61227820 3.1789477 3.50564054 3.3213531
                   Cd1d1      Tex15     Prss35      Esr1     Efcab2      Myzap      Rtp4      Sept3      Nlrc5     Zfp947      Trpc1
Nmu           3.61937291 3.66119462 3.63716950 3.6452166 3.70447020 3.73225185 3.8177538 3.71775212 3.70730065 3.74898656 3.76732421
Fam159b       3.43818716 3.45998343 3.46205583 3.4842244 3.50841948 3.53012454 3.7455479 3.52937397 3.53840951 3.55460761 3.56851541
Them5         4.86223751 5.04673927 4.82897910 4.7271410 5.03963519 5.10787467 4.1553476 4.99032731 4.83246991 5.06035991 5.10880442
Hsd3b2        3.50203253 3.62899157 3.49193441 3.4388986 3.64551889 3.69836750 3.2422279 3.62361785 3.52929768 3.67895012 3.71567438
                    Mxd4     Syngr3     Armcx6    Trim12c      Casz1    Slc17a1      H2-Q5       Irf7     Strip2 5031425E22Rik        Nmi
Nmu           3.76897154 3.75779873 3.74695121 3.77108930 3.76723803 3.77106151 3.77711326 3.81686902 3.83700451    3.84848629 3.84994883
Fam159b       3.57157987 3.57271562 3.58845321 3.58520575 3.63357492 3.59678302 3.61065417 3.68926211 3.63853602    3.65688863 3.65745819
Them5         5.09900562 4.99482912 4.78121038 5.01101951 4.60356581 4.92092752 4.86477513 4.59132487 5.15952831    5.11430138 5.12243724
Hsd3b2        3.71082315 3.64588275 3.52143473 3.66162451 3.43693977 3.61098041 3.58293952 3.45792727 3.77840156    3.75824721 3.76357504
                    Btg2   Arhgap33      Adck2      Dtwd2      Psen2     Rasl12      Csf1r       Pex6      Arl10     Zfp449       Mira
Nmu           3.86744034 4.58063126 4.65814437 4.60079334 4.62787235 4.62237700 4.68510914 4.61565794 4.63725555 4.64001964 4.67408471
Fam159b       3.66641030 4.34555681 4.41013877 4.36477234 4.38647006 4.38241957 4.43531477 4.38016550 4.39898149 4.40402208 4.43973793
Them5         5.20320431 6.06490717 6.24495099 6.09023172 6.16048694 6.14307900 6.28413403 6.09820297 6.14093010 6.12328841 6.13772545
Hsd3b2        3.81809793 4.64778918 4.78664553 4.67161448 4.72424823 4.71175941 4.82154001 4.68310448 4.71746563 4.70874986 4.73293932
                    Eya4       Zxda   BC017158      Folr1      Fcho1    Zdhhc14     Klhl30   Cyp4a12b      Lppr2    Slc35d2       Cbr2
Nmu           4.69051473 4.68321115 4.74483690 4.75266908 4.77563764 4.76334461 4.82959824 4.76513655 4.80221644 4.79723821 4.84828165
Fam159b       4.45381103 4.44947583 4.50147036 4.50883673 4.53237611 4.52269144 4.57909623 4.52658663 4.55747527 4.55412289 4.59800915
Them5         6.17243038 6.14010739 6.27798959 6.28887986 6.30352122 6.26967360 6.41539116 6.25249360 6.33967600 6.32079690 6.42951572
Hsd3b2        4.76026413 4.73859207 4.84549286 4.85532471 4.87431697 4.84939620 4.96302015 4.84052898 4.90720325 4.89417287 4.97967523
                 Ccdc142     Fam98c   Synpo2l        Cpq        Bok     Ubxn11    Rnf144a    Arhgdib     Acss3  Kcnq1ot1    Nkain4
Nmu           4.78184148 4.86317868 5.2942193 4.85814116 4.92689723 4.86632086 4.89593108 4.92438773 4.8482083 5.0890228 5.4178434
Fam159b       4.54513604 4.61811474 5.0022996 4.61509921 4.67739855 4.62807145 4.65260759 4.67684030 4.6170912 4.8220156 5.1248143
Them5         6.25040664 6.39522381 7.2288865 6.37263627 6.49078959 6.33665359 6.40782815 6.47080887 6.2575693 6.7967193 7.3533783
Hsd3b2        4.84722203 4.96704952 5.6438738 4.95192387 5.05081848 4.93552422 4.98945622 5.03831740 4.8828030 5.3000506 5.7704748
                    Hhat  Rab11fip4    Tgfb1i1     Slamf9        Cfp      Acta1    Slc24a3    Kirrel3        Bik    Pitpnm3      Nkd2
Nmu           5.00756977 5.04807053 5.11957195 5.11440782 5.15224764 5.05828232 5.05266595 5.19885423 5.24456522 5.22047694 5.3358788
Fam159b       4.75953160 4.80331689 4.86454396 4.86300261 4.90060182 4.82563884 4.82393812 4.95555754 5.00028408 4.98123767 5.0831683
Them5         6.54771382 6.55293778 6.71061268 6.67231116 6.70781280 6.45179251 6.41170967 6.67131468 6.72051930 6.65286654 6.8799726
Hsd3b2        5.12047580 5.14235680 5.26451901 5.24045985 5.27805434 5.09110307 5.06642214 5.27955147 5.32851952 5.27967470 5.4602788
                     Nog      Acta2    Ccser1   Gm10416     Col7a1       Fbp2    Trabd2b      Upk3b 3110070M22Rik   Gm11423     Fam78b
Nmu           5.27007397 5.24619486 5.4902839 5.3351555 5.47129236 5.65785873 5.50082169 5.72598215    5.75228852 5.4723684 5.78309981
Fam159b       5.02830993 5.00891744 5.2319548 5.1010007 5.21888505 5.38662543 5.25261110 5.45284325    5.48007705 5.2389511 5.51067421
Them5         6.71960127 6.65732341 7.0710678 6.7064371 6.99739517 7.34860329 6.98379483 7.42959267    7.44415221 6.8200437 7.47434057
Hsd3b2        5.34015816 5.29448086 5.6391783 5.3643420 5.58921278 5.87269756 5.59566025 5.94959893    5.96978853 5.4927586 6.00089455
                   Mrvi1     Ttc21a     Rapsn     Icam5  Tbc1d10c       Mak      Nppb   Ccdc116    Atp2a1    Dusp13      Gbp8    Mapk10
Nmu           5.71404886 5.69497638 6.4056892 6.0877699 6.2480414 6.5576248 6.4477428 6.4511044 6.7955177 7.0005041 6.7170580 6.4407350
Fam159b       5.45893126 5.44388470 6.1070140 5.8316148 5.9757583 6.2743336 6.1725145 6.1943888 6.5157108 6.7082039 6.4430384 6.1946227
Them5         7.24062011 7.18458896 8.3304148 7.5885274 7.8983297 8.3013718 8.1132946 7.9257176 8.4852814 8.8161362 8.3501719 7.8102497
Hsd3b2        5.83822716 5.79836576 6.7686204 6.2062356 6.4520331 6.8192763 6.6638346 6.5628620 7.0308012 7.3092273 6.9197995 6.4962459
                    Nrp Serpinb9e     Bmp8b    Mamdc2     Myo7b    Olfr99 4930550C14Rik   Ankrd61     Ephx4    Fcer1g   Mir1191   Mir8091
Nmu           6.4407350 6.6322896 7.1618078 7.0617038 7.1413768 6.9754529     7.1796794 7.3763598 7.3307293 7.5771839 7.4616975 7.3307293
Fam159b       6.1946227 6.3653658 6.8678967 6.8001572 6.8729179 6.7225851     6.9195096 7.1004443 7.0580606 7.2890791 7.1802240 7.0580606
Them5         7.8102497 8.1970608 8.9867009 8.5390903 8.6861745 8.3694191     8.6344390 8.9867009 8.9085201 9.3125433 9.1284449 8.9085201
Hsd3b2        6.4962459 6.7961950 7.4779023 7.1859077 7.3029191 7.0540605     7.2938011 7.5767576 7.5128845 7.8485932 7.6939487 7.5128845
                  Npas3     Smyd1 1810009A15Rik      Il23r  Tnfrsf25      Gulo    Kcnip2     Tssk2     Grin1     Sh3bgr       Lyz2
Nmu           7.3307293 7.4372709     7.5361477  8.4434229 7.8768034 7.7148500 7.7148500 7.8394036 8.1435968  9.1455647  9.0014038
Fam159b       7.0580606 7.1659180     7.2659871  8.1314188 7.5860841 7.4467770 7.4467770 7.5626640 7.8802035  8.8401594  8.7105493
Them5         8.9085201 8.9942144     9.0741435 10.4335516 9.6270537 9.2195445 9.2195445 9.4322200 9.5730413 11.0241397 10.7035238
Hsd3b2        7.5128845 7.6097576     7.6999054  8.8668092 8.1600041 7.8633991 7.8633991 8.0365701 8.2583842  9.5118950  9.2699707
 [ reached getOption("max.print") -- omitted 209 rows ]

> # Make dendrogram data----
> dhc <- as.dendrogram(hclust(d = sampleDists),
+                      horiz = TRUE)

> ddata <- dendro_data(dhc, 
+                      type = "rectangle")

> # Gene orger----
> ddata$labels
      x y         label
1     1 0         Acta1
2     2 0       Slc24a3
3     3 0       Gm10416
4     4 0       Gm11423
5     5 0       Pitpnm3
6     6 0         Acta2
7     7 0           Nog
8     8 0       Kirrel3
9     9 0           Bik
10   10 0         Adck2
11   11 0         Csf1r
12   12 0        Klhl30
13   13 0          Cbr2
14   14 0           Bok
15   15 0       Arhgdib
16   16 0        Fam98c
17   17 0           Cpq
18   18 0       Rnf144a
19   19 0         Lppr2
20   20 0         Fcho1
21   21 0       Slc35d2
22   22 0       Zdhhc14
23   23 0      BC017158
24   24 0         Folr1
25   25 0         Psen2
26   26 0        Rasl12
27   27 0          Eya4
28   28 0          Mira
29   29 0          Zxda
30   30 0      Arhgap33
31   31 0         Dtwd2
32   32 0          Pex6
33   33 0         Arl10
34   34 0        Zfp449
35   35 0         Acss3
36   36 0        Ubxn11
37   37 0      Cyp4a12b
38   38 0       Ccdc142
39   39 0      Kcnq1ot1
40   40 0          Hhat
41   41 0     Rab11fip4
42   42 0       Tgfb1i1
43   43 0        Slamf9
44   44 0           Cfp
45   45 0       Synpo2l
46   46 0        Nkain4
47   47 0          Fbp2
48   48 0         Upk3b
49   49 0 3110070M22Rik
50   50 0        Fam78b
51   51 0         Mrvi1
52   52 0        Ttc21a
53   53 0        Ccser1
54   54 0          Nkd2
55   55 0        Col7a1
56   56 0       Trabd2b
57   57 0         Il23r
58   58 0        Sh3bgr
59   59 0          Lyz2
60   60 0         Grin1
61   61 0         Smyd1
62   62 0 1810009A15Rik
63   63 0       Ankrd61
64   64 0         Npas3
65   65 0         Ephx4
66   66 0       Mir8091
67   67 0         Tssk2
68   68 0          Gulo
69   69 0        Kcnip2
70   70 0      Tnfrsf25
71   71 0        Fcer1g
72   72 0       Mir1191
73   73 0         Rapsn
74   74 0        Dusp13
75   75 0         Bmp8b
76   76 0         Icam5
77   77 0       Ccdc116
78   78 0        Mapk10
79   79 0           Nrp
80   80 0        Olfr99
81   81 0         Myo7b
82   82 0        Mamdc2
83   83 0 4930550C14Rik
84   84 0           Mak
85   85 0        Atp2a1
86   86 0          Gbp8
87   87 0     Serpinb9e
88   88 0      Tbc1d10c
89   89 0          Nppb
90   90 0       Gm19589
91   91 0         Iigp1
92   92 0         Them5
93   93 0        Cd300a
94   94 0        Hsd3b2
95   95 0        Zfp345
96   96 0        Entpd1
97   97 0         Scn1a
98   98 0         Obscn
99   99 0 5330417C22Rik
100 100 0           Nmu
101 101 0       Fam159b
102 102 0        Nckap5
103 103 0         Sox21
104 104 0         Cd59a
105 105 0       Mkln1os
106 106 0         Akap5
107 107 0        Fbxl13
108 108 0           C4b
109 109 0       Fam196a
110 110 0       Snora5c
111 111 0        Kcnj15
112 112 0 4933433G19Rik
113 113 0       Alox12b
114 114 0       Gm15708
115 115 0        Phf11a
116 116 0 D630023F18Rik
117 117 0       Slc35f3
118 118 0        Aknad1
119 119 0         Thbs2
120 120 0        Gm3704
121 121 0       Mab21l3
122 122 0        Slc5a1
123 123 0          Fgl2
124 124 0        Dmrta2
125 125 0         Bambi
126 126 0         Etnk2
127 127 0         Ifi44
128 128 0       Col11a2
129 129 0          Esam
130 130 0        Nipal1
131 131 0         Efr3b
132 132 0         Grem1
133 133 0       Pcdhgc4
134 134 0          Gng2
135 135 0          Rgcc
136 136 0         Epcam
137 137 0          Ism1
138 138 0          Rtp4
139 139 0           Ntm
140 140 0          Pax9
141 141 0         Alcam
142 142 0         Casz1
143 143 0          Irf7
144 144 0         Klra2
145 145 0      Rasgef1a
146 146 0        Inpp5d
147 147 0 A630066F11Rik
148 148 0       Pcdhgb6
149 149 0          Evpl
150 150 0       Fam212a
151 151 0         Adcy7
152 152 0        Celsr1
153 153 0         Mdga1
154 154 0        Dnah7b
155 155 0           Fv1
156 156 0         Ptpn7
157 157 0           Sp6
158 158 0         Tex15
159 159 0         Sept3
160 160 0        Syngr3
161 161 0       Trim12c
162 162 0         Myzap
163 163 0         Trpc1
164 164 0          Mxd4
165 165 0        Efcab2
166 166 0        Zfp947
167 167 0 5031425E22Rik
168 168 0           Nmi
169 169 0        Strip2
170 170 0          Btg2
171 171 0       Slc17a1
172 172 0         Nlrc5
173 173 0         H2-Q5
174 174 0          Esr1
175 175 0        Armcx6
176 176 0       Gm15446
177 177 0          Has3
178 178 0          Oas3
179 179 0         Grhl3
180 180 0        Cyfip2
181 181 0        Prss35
182 182 0          Shc4
183 183 0         Cd1d1
184 184 0        Clstn3
185 185 0       Pip5k1b
186 186 0         Tdrd5
187 187 0        Slc2a4
188 188 0         Kcnj2
189 189 0         Esrp1
190 190 0        Thsd7a
191 191 0          Hey2
192 192 0          Map7
193 193 0          Hic1
194 194 0     Hist1h2ba
195 195 0          Liph
196 196 0        Pkd1l3
197 197 0          Gypc
198 198 0       Pacsin1
199 199 0       Slc6a18
200 200 0        Lrrc15
201 201 0         Ltc4s
202 202 0          Otog
203 203 0         Ptges
204 204 0         Catip
205 205 0        Sprr1b
206 206 0         Ccr10
207 207 0         Itpka
208 208 0       Ptgs2os
209 209 0       Fam110c
210 210 0        Rimbp3
211 211 0          Plb1
212 212 0       Slc6a12
213 213 0         Dbil5

> # Segment data----
> dtp1 <- segment(ddata)

> # Hitmap data----
> dtp2 <- ll

> dtp2$Gene <- factor(dtp2$Gene,
+                     levels = ddata$labels$label)

> offset.size <- 10

> p1 <- ggplot(data = dtp2) +
+   coord_polar("y",
+               start = 0,
+               direction = -1) +
+   geom_tile(aes(x =  as.numeric(Comp .... [TRUNCATED] 

> p1

> tiff(filename = "tmp/MES13_RNA_DEGseq_TIIA-HG-LG_hitmap_with_phylo.tiff",
+      height = 12,
+      width = 12,
+      units = 'in',
+      res = 1 .... [TRUNCATED] 

> plot(p1)

> graphics.off()

> # CHECK:
> dt1.fpm[GeneNames == "Iigp1"]
   GeneNames    LG_fpm   HG_fpm TIIA_5uM_fpm
1:     Iigp1 0.9534629 2.265238    0.1028917

> dtp2[Gene == "Iigp1"]
    Gene Comparison Gene Expression Diff
1: Iigp1      HG-LG             1.248414
2: Iigp1    TIIA-HG            -4.460464

> # Correct
> 
> sessionInfo()
R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggdendro_0.1-20             knitr_1.22                  DEGseq_1.38.0               qvalue_2.16.0              
 [5] DESeq2_1.24.0               SummarizedExperiment_1.14.0 DelayedArray_0.10.0         BiocParallel_1.18.0        
 [9] matrixStats_0.54.0          Biobase_2.44.0              GenomicRanges_1.36.0        GenomeInfoDb_1.20.0        
[13] IRanges_2.18.0              S4Vectors_0.22.0            BiocGenerics_0.30.0         ggplot2_3.1.1              
[17] data.table_1.12.2          

loaded via a namespace (and not attached):
 [1] bit64_0.9-7            splines_3.6.0          Formula_1.2-3          latticeExtra_0.6-28    blob_1.1.1            
 [6] GenomeInfoDbData_1.2.1 pillar_1.3.1           RSQLite_2.1.1          backports_1.1.4        lattice_0.20-38       
[11] digest_0.6.18          RColorBrewer_1.1-2     XVector_0.24.0         checkmate_1.9.3        colorspace_1.4-1      
[16] htmltools_0.3.6        Matrix_1.2-17          plyr_1.8.4             XML_3.98-1.19          pkgconfig_2.0.2       
[21] genefilter_1.66.0      zlibbioc_1.30.0        xtable_1.8-4           scales_1.0.0           htmlTable_1.13.1      
[26] tibble_2.1.1           annotate_1.62.0        withr_2.1.2            nnet_7.3-12            lazyeval_0.2.2        
[31] survival_2.44-1.1      magrittr_1.5           crayon_1.3.4           memoise_1.1.0          MASS_7.3-51.4         
[36] foreign_0.8-71         tools_3.6.0            stringr_1.4.0          munsell_0.5.0          locfit_1.5-9.1        
[41] cluster_2.0.9          AnnotationDbi_1.46.0   compiler_3.6.0         rlang_0.3.4            grid_3.6.0            
[46] RCurl_1.95-4.12        rstudioapi_0.10        htmlwidgets_1.3        labeling_0.3           bitops_1.0-6          
[51] base64enc_0.1-3        gtable_0.3.0           DBI_1.0.0              reshape2_1.4.3         gridExtra_2.3         
[56] bit_1.1-14             Hmisc_4.2-0            stringi_1.4.3          Rcpp_1.0.1             geneplotter_1.62.0    
[61] rpart_4.1-15           acepack_1.4.1          xfun_0.6              

> sink()
