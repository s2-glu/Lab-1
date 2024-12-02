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
10. These commands are from #3 which links a table. When you have finished, answer the following questions:
