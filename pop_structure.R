require(RColorBrewer)
palette(brewer.pal(12, "Set3"))

plot_freq_vs_HWE<-function(file,title="",                                        
    cols=c("#ffb000", "#89bff7", "#004400")){                                    
                                                                                 
    plink_hwe <- read.table(file, as.is=TRUE, header=T)                          
                                                                                 
    counts <- sapply(plink_hwe$GENO,function(x){as.numeric(strsplit(x,"/")[[1]])})
        counts<-t(counts)                                                        
        tot_counts<-rowSums(counts)                                              
        geno_freq<-counts / tot_counts                                           
        allele_freq<-(geno_freq[,1]+.5 *geno_freq[,2])                           
                                                                                 
        ss_allele<-c(allele_freq,1-allele_freq)                                  
        ss_geno<-rbind(geno_freq,geno_freq[,c(3,2,1)])                           
                                                                                 
    col_alpha <- adjustcolor(cols, alpha.f=0.05)
    allele_jitter <- jitter(ss_allele, 1.5)
    ss_geno <- jitter(ss_geno, 1.5)
    ss_geno <- pmax(pmin(ss_geno, 1), 0)
    allele_jitter <- pmax(pmin(allele_jitter, 1), 0)
    plot(allele_jitter,ss_geno[,1],xlim=c(0,1),ylim=c(0,1.5),col=col_alpha[1],         
        xlab="allele frequency",ylab="genotype frequency",main=title, pch=16)            
    points(allele_jitter,ss_geno[,3],xlim=c(0,1),ylim=c(0,1),col=col_alpha[2], pch=16)       
    points(allele_jitter,ss_geno[,2],xlim=c(0,1),ylim=c(0,1),col=col_alpha[3], pch=16)       
                                                                                 
    smooth=1/5                                                                   
    lines(lowess(ss_geno[,1]~ss_allele,f = smooth),col="black", lwd=4)           
    lines(lowess(ss_geno[,3]~ss_allele,f = smooth),col="black", lwd=4)           
    lines(lowess(ss_geno[,2]~ss_allele,f = smooth),col="black", lwd=4)           
                                                                                 
        x=1:1000/1000                                                            
        lines(x,x^2,lty=2, lwd=4)                                                       
        lines(x,2*x*(1-x),lty=2, lwd=4)                                                 
        lines(x,(1-x)^2,lty=2, lwd=4)                                                   
        legend(x=0.0,y=1.5,col=c(cols,rep("black",2)),                             
            legend=c("Homozygote AA","Homozygote aa","Heterozygote Aa",          
                "Mean","Hardy Weinberg Expectation"),pch=c(rep(16,3),            
                    rep(NA,2)),lty=c(rep(NA,3),1,2))                             
}

load_pca <- function(pc_file, fam_file, name_file){
    require(dplyr)
    pcs <- read.table(pc_file)
    names(pcs) <- sprintf("PC%s",  1:ncol(pcs))
    fam <- read.table(fam_file, as.is=T)[,1:2]
    names(fam) <- c('IID', 'FID')
    name <- read.table(name_file, as.is=T)
    names(name) <- c('IID', 'abbrev', 'order')
    return(fam %>% left_join(name) %>% cbind(pcs))
}

plot_pca <- function(pc_data, x_axis='PC1', y_axis='PC2', ...){
    plot(pc_data[,x_axis], pc_data[,y_axis], type='n',...)
    text(pc_data[,x_axis], pc_data[,y_axis], label=pc_data$abbrev)
}

load_admixture_ll <- function(K=2:12, reps=1:10, template_names="admixture/%s/D_magna3_%s.like"){
    liks <- sapply(K, function(k) sapply(reps, function(i) read.table(sprintf(template_names, i, k))[,2]))
    colnames(liks) <- sprintf("K%s", K)

    as.data.frame(liks)
}

sort_admixture_Q <- function(adm, n){
    adm <- cbind(n, adm)
    adm %>% arrange(order)
}

load_admixture_Q_single_K <- function(K=2, reps=1:10, template_names="admixture/%s/D_magna3.%s.Q", name_file='data/D_magna.names'){
    n <- read.table(name_file)
    names(n) <- c('IID', 'abbrev', 'order')
    lapply(reps, function(i){ adm <- read.table(sprintf(template_names, i, K))
                adm <- sort_admixture_Q(adm, n)
            })
}

load_admixture_Q_single_rep <- function(K=2:12, template_names="admixture/D_magna3.%s.Q", name_file='data/D_magna.names'){
    n <- read.table(name_file)
    names(n) <- c('IID', 'abbrev', 'order')
    lapply(K, function(i){ adm <- read.table(sprintf(template_names, i))
                adm <- sort_admixture_Q(adm, n)
            }
           )
}
load_admixture_Q_best <- function(lik, K=2:5, reps=1:20, template_names="admixture/%s/D_magna3.%s.Q", name_file='data/D_magna.names'){
    n <- read.table(name_file)
    names(n) <- c('IID', 'abbrev', 'order')

    best_rep_per_K <- apply(lik, 2, which.max)

    mapply( function(k, i){
	    adm <- read.table(sprintf(template_names, i, k))
	    adm <- sort_admixture_Q(adm, n)
	   }, K, best_rep_per_K)
}

plot_admixture_single <- function(df, ...){
	barplot(t(df[,-1:-3]), col=1:12,
		names.arg=df$abbrev, las=2, ...)
}

