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
Explanation: 
1. Use a  startingproteinprotein from Homo sapiens (human) as the query sequence. It can be downloaded using the following command.
2. Perform a blast search using the query protein.
3. Look at the output in globins.blastp.typical.out using less
4. Create a more detailed and easier-to-process output of the same analysis. The -outfmt flag specifies a particular output format that will be useful for our analysis. 
5. Look at the output in startingprotein.blastp.detail.out using the less -S command.
6. Use this command to avoid counting by hand.
7. Use this command you to filter our output file to satisfy this requirement.
8. Use wc command to count the total number of hits in the BLAST results after the filter.
9. Use grep to know how may paralogs are found in each species.
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
Explanation:
1. Can grab the sequences we want (their names are in globins.blastp.detail.filtered.out) from allprotein.fas using seqkit. This is the seqkit command to obtain the sequences that are in the BLAST output file.
2. We want to align all of the sequences with each other along their entire length. The resulting alignment will be a hypothesis about which positions are homologous to each other among the globin proteins, and which positions contain insertions or deletions, with respect to the other sequences. Use command to make a multiple sequence alignment using muscle.
3. Use the R package msa and a script that we have provided for you to print your alignment to a large pdf file that you can inspect with ease. Run this Rscript for printing the alignment.
4. Calculate the width (length) of the alignment.
5. Calculate the length of the alignment after removing any column with gaps.
6. Calculate the length of the alignment after removing invariant (completely conserved) positions.
7. Calculate the average percent identity using t_coffee.
8. Repeat calculating the average percent identity using alignbuddy.
   
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
Explanation:
1. Create a directory in lab 5 for the startingprotein tree.
2. Go to that directory.
3. Run the following command, which will both remove any sequence that contains a duplicate label tag, and put a copy in the lab05 directory.
4. Use IQ-TREE to find the maximum likehood tree estimate. First, it will calculate the optimal amino acid substitution model and amino acid frequencies. Then, it will perform a tree search, estimating branch lengths as it goes.
5. The .iqtree file includes an ASCII graphics (text graphics) version of the tree. It can also display it by reading the .treefile (which is newick formatted) into the nw_display program.
6. look at it unrooted with a graphical display. nw_display can't do this, use our R script instead. Because there are so many genes,  the size of the text labels is smaller (0.4), and set the label lengths to 15.
7. Use the software gotree to reroot the tree by using a type of rooting called midpoint.
8. Look at the rooted tree at the command line, but making it graphic instead.
9. Convert this svg image to a pdf.
10. The tree shown by default in nw_display (and most other programs) is a phylogram. This means that the lengths of each branch are proportional to the number of substitutions that have accumulated in the sequence along that branch. This command witching the view to a cladogram
# Lab-6
Order:
1. mkdir ~/lab06-s2-glu/startingprotein
2. cp ~/lab05-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.al.mid.treefile
3. java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-s2-glu/species.tre -g ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-s2-glu/startingprotein/
4. nw_display ~/lab05-s2-glu/species.tre
5. grep NOTUNG-SPECIES-TREE ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.ntg | sed -e "s/^[&&NOTUNG-SPECIES-TREE//" -e "s/]/;/" | nw_display -
6. gotree prune -i ~/lab06-s2-glu/startingprotein/startingprotein.homologs f.al.mid.treefile -o ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H. sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta
7. python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.ntg --include.species
8. thirdkind -Iie -D 40 -f ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.ntg.xml -o ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.svg
9. convert -density 150 ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.svg ~/lab06-s2-glu/startingprotein/startingprotein.homologsf.pruned.treefile.rec.pdf
Explanation:
1. Create a new folder in lab06 for startingprotein.
2. Changing the name mygenefamily to match startingprotein. Then, make a copy of this gene tree from lab 5 startingprotein folder into your lab 6 startingprotein folder.
3. Command to perform the reconciliation.
4. Note the provided names for internal lineges.
5. See the node names that notung used/assigned to these internal nodes with this command.
6. Prune the tree to see important lineages.
7. Generate a RecPhyloXML object by running the python script.
8. Use software that allows us to view the gene tree reconciliation within the species tree.
9. Convert this to a pdf.
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
Explanation:
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
