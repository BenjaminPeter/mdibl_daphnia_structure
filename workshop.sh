#download
wget -O data/D_magna_RADseq.vcf.gz http://datadryad.org/bitstream/handle/10255/dryad.92017/D_magna_RADseq.vcf.gz 

#preprocessing
plink --vcf data/D_magna_RADseq.vcf.gz --make-bed --out data/D_magna --allow-extra-chr --geno 0.0  --remove data/to-remove.txt 
ln -f data/D_magna.bed data/D_magna2.bed 
ln -f data/D_magna.fam data/D_magna2.fam 
awk '{print 1, $2, $3, $4, $5, $6}' data/D_magna.bim > data/D_magna2.bim                                                           
plink --bfile data/D_magna2 --recode A --out data/D_magna2                                                                         
cut -f1,7- -d" " data/D_magna2.raw | tail -n+2 > data/D_magna.txt                                                                  

# hwe test
plink --bfile data/D_magna --hardy   --out

#filter
plink --bfile data/D_magna2 --hwe 0.001 --make-bed --out data/D_magna3
plink --bfile data/D_magna3 --hardy --out data/D_magna3

mkdir -p admixture pca

#run flashpca
flashpca --bfile data/D_magna3 --outpc pca/D_magna3_pc.txt --outpve pca/D_magna3_pve.txt --outload pca/D_magna3_load.txt   --ndim 22

#run admixture
for i in `seq 2 8`; do
	for j in `seq 1 20`; do
		mkdir -p admixture/$j/
		ln -f data/D_magna3.bed admixture/$j
		ln -f data/D_magna3.bim admixture/$j
		ln -f data/D_magna3.fam admixture/$j
		cd admixture/$j
		echo "K=$i/ run $j"
		admixture D_magna3.bed $i --seed=$j | grep ^Log > D_magna3_$i.like
		cd ../..
	done
done

