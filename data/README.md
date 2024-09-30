# all_meta_corrected_host.csv

Metadata file for ssRNA viruses included in the analysis (see results section), with the column names defined as follows:

Corrected_host: Corrected virus host category

Strain: Virus strain name

Accession: Virus accession number

Date: Virus collection date

Region: Virus collection region

Host: Original recorded host for the virus

Species: Virus species name

Family: Virus family name

Order: Virus order name

Year: Virus collection year

 

# all_mpxv_sars_phyvirus_sample.txt

The 192 mutation types for ssRNA viruses included in the analysis (see results section). For example, if the mutation type is C>U, with A at the 5’ end and G at the 3’ end, it is represented as CUAG.

 

# all_signature_sample_csv.txt

The 192 mutation types for ssRNA viruses included in the analysis (see results section), aggregated by time and host category.

 

# all_phyvirus_human_compos.txt

The possible mutation types for ssRNA viruses and the human genome. 

 

# context_all_signature_sample_csv.txt

The 192 mutation types for ssRNA viruses, adjusted according to the human genome (see methods).

 

# context_96_all_signature_sample.txt

The mutation spectrum from context_all_signature_sample_csv.txt, folded into 96 mutation types.

 

# SBS96_De-Novo_Signatures.txt

Nine de novo signatures extracted from context_96_all_signature_sample.txt using SigProfilerExtractor.

 

# SBS96_De-Novo_Activities_refit.txt

The number of mutations in each ssRNA virus mutation spectrum explained by the 9 de novo signatures.

 

# CH192_De-Novo_Signatures_sars.txt

Two de novo signatures extracted from SARS-CoV-2 using SigProfilerExtractor.

 

# CH192_De-Novo_Activities_refit_sars.txt

The number of mutations in each SARS-CoV-2 mutation spectrum explained by the 2 de novo signatures.

 

# sars_coding_mutation.txt

The mutation file for the SARS-CoV-2 coding region, with the column names defined as follows:

accession: Virus accession number

prot_pos: Amino acid position of the mutation

real_pos: nucleotide position of the mutation

sub_type: Mutation type (synonymous or non-synonymous)

mutation_type: nucleotide mutation type

Sub: amino acid substitution type

 

# sars_exp_site.txt

The expected number of synonymous and non-synonymous mutation sites in SARS-CoV-2.

 

# CH192_De-Novo_Signatures_h3n2.txt

Three de novo signatures extracted from H3N2 using SigProfilerExtractor.

 

# CH192_De-Novo_Activities_refit_h3n2.txt

The number of mutations in each H3N2 virus mutation spectrum explained by the 3 de novo signatures.

 

# h3n2_coding_mutation.txt

The mutation file for the H3N2 coding region.

 

# H3n2_exp_site.txt

The expected number of synonymous and non-synonymous mutation sites in H3N2.

 

# all_pairwise_result1.csv

Cosine similarity of mutation spectra between iSNVs and SNPs in SARS-CoV-2.

 

# syn_pairwise_result2.csv

Cosine similarity of synonymous mutation spectra between iSNVs and SNPs in SARS-CoV-2.

 

# nonsyn_pairwise_result3.csv

Cosine similarity of non-synonymous mutation spectra f between iSNVs and SNPs in SARS-CoV-2.

 

# isnv_pair_ruan_andreas.csv

iSNV transmission file for SARS-CoV-2 , with the column names defined as follows:

donor: Donor ID

recipient: Recipient ID

donor_alt_frequency: Frequency of the iSNV in the donor

recipient_alt_frequency: Frequency of the iSNV in the recipient

func: functional changes resulted from the iSNV 

ref_pos_alt: Type of iSNV mutation

recipient_count: Number of recipients for the donor