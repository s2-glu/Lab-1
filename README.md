# Lab-3
Order:
1. ncbi-acc-download -F fasta -m protein "NP_001914.3"
2. blastp -db ../allprotein.fas -query NP_001914.3.fa -outfmt 0 -max_hsps 1 -out globins.blastp.typica l.out
3. less startingprotein.blastp.typical.out
4. blastp -db ../allprotein.fas -query NP_001914.3.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out startingprotein.blastp.detail.out
5. less -S startingprotein.blastp.detail.out
6. grep -c H.sapiens startingprotein.blastp.detail.out
7. awk '{if ($6< 1e-30)print $1 }' startingprotein.blastp.detail.out > startingprotein.blastp.detail.filtered.out
8. wc -l startingprotein.blastp.detail.filtered.out
9. grep -o -E "^[A-Z].[a-z]+" startingprotein.blastp.detail.filtered.out | sort | uniq -c

# Lab-4
Order:
1. seqkit grep --pattern-file ~/lab03-s2-glu/startingprotein/startingprotein.blastp.detail.filtered.out ~/lab03-s2-glu/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-s2-glu/startingprotein/startingprotein.homologs.fas
2. muscle -align ~/lab04-s2-glu/startingprotein/startingprotein.homologs.fas -output ~/lab04-s2-g lu/startingprotein/startingprotein.homologs.al.fas -alv -kli ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas | less -RS -alv -kli --majority ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas | less -RS
3. Rscript --vanilla ~/lab04-s2-glu/plotMSA.R ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas
4. alignbuddy -al ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas
5. alignbuddy -trm all ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas | alignbuddy -al
6. alignbuddy -dinv 'ambig' ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas | alignbuddy -al
7. t_coffee -other_pg seq_reformat -in ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas -output sim
8. alignbuddy -pi ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas | awk ' (NR>2) { for (i=2;i<=NF ;i++){ sum+=$i;num++} } END{ print(100*sum/num) } '

# Lab-5
Order:
1. mkdir ~/lab05-s2-glu/startingprotein
2. cd ~/lab05-s2-glu/startingprotein
3. sed 's/ /_/g' ~/lab04-s2-glu/startingprotein/startingprotein.homologs.al.fas | seqkit grep -v -r -p "dupelabel" > ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.fas
4. iqtree -s ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.fas -bb 1000 -nt 2
5. nw_display ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.fas.treefile
6. Rscript --vanilla ~/lab05-s2-glu/plotUnrooted.R ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.fas.treefile ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.fas.treefile.pdf 0.4 15
7. gotree reroot midpoint -i ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.fas.treefile -o ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile
8. nw_order -c n ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile | nw_display - nw_order -c n ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s > ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile.svg -
9. convert ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile.svg ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile.pdf
10. nw_order -c n ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.midCl.treefile.svg -

# Lab-6
Order:
1. mkdir ~/lab06-$MYGIT/mygenefamily
2. cp ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile
3. java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-s2-glu/species.tre -g ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-s2-glu/startingprotein/
4. nw_display ~/lab05-s2-glu/species.tre
5. grep NOTUNG-SPECIES-TREE ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.ntg | sed -e "s/^[&&NOTUNG-SPECIES-TREE//" -e "s/]/;/" | nw_display -
6. gotree prune -i ~/lab06-s2-glu/startingprotein/startingprotein.homologs f.al.mid.treefile -o ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H. sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta
7. python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.ntg --include.species
8. thirdkind -Iie -D 40 -f ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.ntg.xml -o ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.svg
9. convert -density 150 ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.svg ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.pdf

# Lab-8
Order:
1. mkdir ~/lab08-s2-glu/startingprotein && cd ~/lab08-s2-glu/startingprotein
2. sed 's/*//' ~/lab04-s2-glu/startingprotein/startingprotein.homologs.fas > ~/lab08-s2-glu/startingprotein/startingprotein.homologs.fas
3. rpsblast -query ~/lab08-s2-glu/startingprotein/startingprotein.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-s2-glu/startingprotein/startingprotein.rps-blast.out -outfmt "6 qseqid qlen  qstart qend evalue stitle" -evalue .0000000001
4. cp ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.outgroupbeta.treefile ~/lab08-s2-glu/startingprotein
5. Rscript --vanilla ~/lab08-s2-glu/plotTreeAndDomains.r ~/lab08-s2-glu/startingprotein/startingprotein.homologsf.outgroupbeta.treefile ~/lab08-s2-glu/startingprotein/startingprotein.rps-. blast.out ~/lab08-s2-glu/startingprotein/startingprotein.tree.rps.pdf
6. mlr --inidx --ifs "\t" --opprint cat ~/lab08-s2-glu/startingprotein/startingprotein.rps-blast.out | tail -n +2 | less -S
7. cut -f 1 ~/lab08-s2-glu/startingprotein/startingprotein.rps-blast.out | sort | uniq -c
8. cut -f 6 ~/lab08-s2-glu/startingprotein/startingprotein.rps-blast.out | sort | uniq -c
9. awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-s2-glu/startingprotein/startingprotein.rps-blast.out | sort -k2nr
10. cut -f 1,5 -d $'\t' ~/lab08-s2-glu/startingprotein/startingprotein.rps-blast.out

1. First, make a directory for the startingprotein sequences, and change into that directory.
2. Make a copy of our raw unaligned sequence, removing the asterisk (stop codon) in the process. To do this,use sed's substitute command to substitute any instance of an asterisk with nothing. Direct the ouput to startingprteoin folder in lab8.
3. Command to run rps blast, which shows the name of the species, the gene or protein related to their subunit, numbers that can represent sequence lengths and alignment scores, E-value from a sequence alignment, and specific protein family in the Pfam database.
4. Use our final gene tree from lab 5. Copy it over here, as well using cp.
5. Run the script that shows the predicted domains annotated by the tips of your tree, while also seeing a legend with the names of all the PFAM domains.
6. The tab delimited annotations in startingprotein.rps-blast.out. Use mlr to look at it. It produces the query sequence ID, query sequence length, position of the start of the pfam domain in the query sequence, position of the end of the pfam domain in the query sequence, the e-value of the domain prediction, the title of the predicted domain.
7. Command to do this at the command line to see protein annotation.
8. Command to see what Pfam domain annotation is most commonly found.
9. Command to see which protein has the longest annotated protein domain.
10. Command to pull out just the e-values to see which protein has a domain annotation with the best e-value.

