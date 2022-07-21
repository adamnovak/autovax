# Autovax: automatic SARS-CoV-2 Vaccine Candidate Generator

Remember when the first mRNA COVID-19 vaccines were designed in days? Now, with **Autovax**, YOU can design an up-to-date mRNA vaccine for SARS-CoV-2 in *minutes*. Do you want a vaccine candidate targeting Omicron? BA.2? BA.5? BA.2.75? [You can do anything with **Autovax**](https://www.zombo.com/).

## How?

[*Omicron-specific mRNA vaccination alone and as a heterologous booster against SARS-CoV-2*, by Fang *et al.*](https://www.nature.com/articles/s41467-022-30878-4) describes how some researchers at Yale designed, manufactured, and tested their own lipid nanoparticle mRNA (mRNA-LNP) vaccine candidate in mice.

Autovax attempts to automate my understanding of their mRNA design process (replace one block of four amino acids, and change 6 other ones to Proline). It lets you apply it to the spike protein from any full SARS-CoV-2 sequence that has a good copy, so you can design a vaccine for the latest and greatest(?) SARS-CoV-2 variants. It uses BioPython, [Codon Harmony](https://codon-harmony.readthedocs.io/en/latest/readme.html) You can use it to design your own *personal* SARS-CoV-2 vaccine candidate in a few minutes.

## How do I install this?

```
git clone https://github.com/adamnovak/autovax
cd autovax
virtualenv --python python3.9 venv
. venv/bin/activate
pip install -r requirements.txt
```

## Where do I get sequences? I don't want to sell my soul to GISAID!

You can get [the latest SARS-CoV-2 sequences from NCBI](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049). These sequences are not subject to GISAID's terms, and Autovax doesn't use any GISAID sequences itself, so no GISAID data need be harmed in the making of your vaccine candidate.

To get a sequence:

1. Find a row in the table where "Release Date" is recent and "Pangolin" is the variant you fear most accutely.
2. Click on the "Accession". A panel should pop up; you may need to disable your adblocker to see it.
3. Click on the accession again in the panel to go to NCBI's page for the sample.
4. Click on "FASTA" a couple lines below the heading.
5. On the resulting page, copy all the preformatted text starting with `>` and save it to a text file with a `.fa` extension.

## How can I design my own SARS-CoV-2 vaccine candidate for cheaper than the leading brand?

Run the main script against your downloaded .fa file, with the virtualenv you made activated.

```
./autovax.py your-favorite-virus.fa
```

If it complains with one sequence, try another; not all sequences in NCBI have sufficient accuracy in the spike protein to generate a vaccine candidate.

Make sure to run it from the current directory; it needs to run there to find the canonical spike sequence file.

It should print a report like this (which I got by running it on the included canonical spike (`./autovax.py spike.fa`):

```
RNA DETAILS

RNA Name: Vaccine candidate derived from wuhCor1_ncbiGeneBGP
Modifications: 'N1-Methylpseudouridin' [sic]
Type of RNA: Capped mRNA
Template Source: TriLink Supplied
ORF Sequence: ATGTTCGTGTTCCTCGTGCTGCTGCCTCTGGTCTCGTCTCAATGTGTTAACTTGACTACTCGGACCCAGCTGCCACCTGCATACACGAACAGCTTCACCCGAGGAGTTTACTACCCCGACAAAGTATTTCGATCCTCTGTGTTGCATTCCACTCAAGATCTATTCTTGCCATTTTTCTCAAATGTCACTTGGTTTCACGCCATACACGTGTCCGGAACTAACGGCACAAAGCGCTTTGACAATCCAGTGCTGCCCTTCAATGACGGCGTGTACTTCGCCTCCACCGAGAAGAGCAACATTATTCGCGGTTGGATATTCGGGACAACCCTGGACTCCAAGACACAGAGCCTGTTGATTGTGAATAACGCAACAAACGTGGTCATTAAGGTTTGCGAGTTTCAGTTTTGCAATGATCCTTTCCTGGGAGTGTACTACCACAAGAACAACAAGTCTTGGATGGAATCCGAGTTTCGCGTTTATTCTTCTGCCAACAACTGCACGTTTGAATATGTCTCTCAACCCTTCCTGATGGATCTTGAAGGCAAGCAAGGGAACTTCAAGAACCTGAGGGAGTTCGTCTTTAAGAACATCGATGGCTATTTCAAGATTTACTCCAAGCATACGCCAATAAACCTGGTCCGCGATCTACCGCAAGGGTTTTCCGCTCTCGAACCACTTGTGGATCTTCCTATTGGCATTAACATTACAAGATTCCAGACCCTGCTCGCCCTCCATAGATCTTATCTCACACCTGGAGATTCCTCATCCGGCTGGACTGCGGGGGCTGCGGCATACTATGTCGGGTATCTGCAACCACGCACATTCCTCCTGAAGTACAACGAAAATGGAACCATCACTGATGCCGTGGACTGCGCACTCGATCCGCTGAGCGAGACCAAGTGCACTCTCAAGAGTTTCACGGTAGAGAAGGGAATATACCAGACTTCCAATTTCCGGGTTCAGCCAACCGAATCGATAGTCCGTTTCCCAAACATCACTAACTTGTGTCCTTTCGGAGAAGTGTTCAATGCTACTCGCTTCGCCTCTGTTTATGCATGGAATAGAAAGAGGATCAGCAACTGCGTGGCAGATTATAGCGTTCTGTATAACTCTGCATCCTTCTCCACATTCAAGTGTTACGGTGTCTCTCCCACAAAGCTGAACGATCTCTGTTTCACTAACGTCTACGCCGATAGCTTTGTGATCAGGGGCGACGAGGTTAGACAAATCGCACCCGGTCAAACAGGCAAGATAGCAGACTACAACTACAAGCTCCCTGACGATTTCACGGGATGCGTCATAGCCTGGAACAGCAATAACCTGGACTCTAAAGTCGGGGGGAATTATAACTATCTCTACCGTCTTTTCCGAAAATCAAATCTCAAACCTTTTGAACGCGACATCAGCACTGAGATTTATCAAGCGGGATCCACCCCATGCAACGGTGTCGAGGGATTCAATTGCTACTTTCCACTGCAAAGCTATGGGTTCCAACCTACTAATGGGGTCGGGTACCAACCGTATAGAGTAGTCGTCCTGAGCTTCGAGCTGCTCCACGCTCCTGCTACCGTTTGCGGTCCCAAGAAGAGTACCAACCTCGTGAAGAACAAATGTGTCAATTTCAATTTTAACGGGTTGACAGGTACCGGCGTTCTGACTGAGAGTAATAAGAAGTTCCTGCCTTTCCAACAATTCGGGCGTGACATTGCTGACACTACCGACGCTGTTAGGGACCCCCAAACGCTGGAGATTCTCGACATCACTCCATGTAGCTTCGGAGGCGTGTCCGTGATCACACCTGGGACCAACACCAGTAACCAGGTCGCTGTCCTGTATCAAGATGTAAATTGTACTGAAGTTCCTGTTGCCATCCACGCGGACCAGCTTACACCGACTTGGAGGGTCTACTCCACCGGCTCCAACGTCTTCCAAACTAGGGCAGGCTGCCTGATTGGGGCCGAGCATGTCAACAATTCCTACGAATGCGATATTCCCATCGGAGCCGGAATATGCGCATCTTACCAAACTCAGACTAACAGTCCTGGTTCTGCCTCAAGTGTCGCTTCGCAATCCATTATCGCATACACCATGTCCCTGGGAGCCGAGAACAGTGTTGCCTATTCCAACAATAGCATTGCTATCCCTACTAACTTTACCATCTCCGTGACTACCGAAATCCTTCCCGTGTCGATGACCAAGACATCTGTGGATTGCACTATGTACATTTGCGGCGATAGCACCGAGTGCTCCAACCTTCTGCTGCAATATGGCTCTTTCTGCACGCAATTAAATAGAGCGTTGACTGGCATCGCCGTGGAGCAGGACAAGAACACTCAGGAAGTCTTTGCTCAAGTGAAACAAATCTACAAGACCCCACCAATCAAGGACTTTGGTGGTTTTAACTTCTCTCAAATCCTACCCGATCCCTCTAAACCTAGCAAGCGCAGCCCCATCGAGGACCTTCTGTTCAATAAGGTGACACTGGCCGACGCTGGCTTCATTAAACAATATGGAGACTGCCTCGGCGACATTGCCGCTAGAGACCTGATCTGTGCACAAAAGTTCAACGGGCTCACCGTTCTCCCTCCACTCCTGACCGACGAAATGATTGCCCAATATACGTCCGCTCTTCTCGCCGGGACCATCACCTCTGGATGGACATTCGGGGCTGGACCTGCTCTGCAAATCCCATTTCCGATGCAAATGGCATATCGGTTCAACGGCATTGGCGTGACACAGAATGTCCTTTATGAGAATCAAAAGCTGATCGCCAATCAGTTCAATTCTGCGATCGGGAAAATTCAGGACAGTCTTTCTAGTACACCATCTGCTCTGGGGAAACTGCAAGACGTCGTGAACCAAAATGCACAGGCCCTGAACACTCTCGTGAAGCAGCTCTCTTCCAACTTCGGTGCAATTAGCTCCGTTTTGAATGATATTCTCTCTAGGCTGGACCCACCTGAGGCTGAAGTGCAGATTGATCGGCTGATCACTGGCCGCCTGCAATCCTTGCAAACCTACGTGACACAACAACTCATTAGAGCTGCTGAGATTAGGGCTTCTGCTAACCTTGCGGCCACCAAGATGTCTGAATGCGTTCTGGGGCAAAGCAAAAGGGTCGATTTCTGCGGAAAGGGATACCACCTTATGTCATTCCCCCAATCGGCGCCTCATGGCGTGGTTTTCCTTCACGTAACCTATGTGCCGGCGCAAGAAAAAAACTTCACTACCGCTCCCGCTATCTGTCACGACGGCAAGGCTCATTTCCCCCGGGAGGGAGTGTTTGTGTCCAACGGCACCCACTGGTTCGTCACTCAACGAAACTTCTACGAACCTCAGATCATAACCACGGATAACACTTTTGTGTCGGGAAACTGTGACGTTGTCATCGGGATCGTGAACAATACCGTCTATGACCCACTCCAACCCGAACTGGATAGCTTTAAAGAGGAGCTTGACAAGTACTTTAAGAATCACACATCACCTGATGTGGACCTGGGTGACATTTCTGGAATCAACGCAAGCGTTGTGAACATCCAAAAGGAAATAGACCGGCTGAACGAGGTGGCCAAGAACCTCAATGAGAGCCTCATCGATCTGCAAGAGCTGGGCAAATATGAGCAGTATATTAAATGGCCCTGGTACATCTGGCTGGGCTTTATCGCCGGCCTTATCGCAATAGTTATGGTGACAATTATGCTGTGCTGCATGACATCATGCTGTTCTTGCCTGAAGGGGTGCTGCTCATGTGGCAGCTGTTGCAAGTTCGATGAGGACGATTCTGAACCTGTTCTCAAAGGAGTCAAATTACACTACACATAG

MANUFACTURING PROCEDURE

1. Enter above details into <https://www.trilinkbiotech.com/mrna.html> and order.
2. Assemble into lipid nanoparticles with NanoAssemblr <https://www.precisionnanosystems.com/platform-technologies/product-comparison/ignite>.
3. Buffer exchange or something, idk.
4. Inject into your (least?) favorite mice.

For more information, see <https://www.nature.com/articles/s41467-022-30878-4>

Thank you for using Autovax 🌈
```

## What should I do with my vaccine candidates?

No immunologists or actual wet lab biologists were consulted in the development of Autovax, so the candidates and instructions it produces might not work at all, and might kill or injure any mice, ferrets, humans, or other animals who are subjected to them. They *definitely* have not been through any sort of clinical trial, nor has this code, and I am probably required by law to tell you that **nothing here has been evaluated by the FDA** and **nothing here is intended to cure, treat, or prevent any disease**. Go take an actually-tested vaccine!

I have not contacted the authors of the paper I am working from to get all of the no-doubt essential information (such as their actual construct sequences) that I could not get from the paper and its supplement. I borrowed the "Hexapro" substitutions, which seem to have a brand name for six Prolines, from a web site linked from a paper that seemed to purposefully leave out the actual substitutions to make; they may be patented or something to the extent that you can patent six Prolines. I have no idea how you get an mRNA into an LNP, or if the mRNA provider that Autovax suggests can actually produce the mRNA you would need. If you want real help replicating [*Omicron-specific mRNA vaccination alone and as a heterologous booster against SARS-CoV-2*, by Fang *et al.*](https://www.nature.com/articles/s41467-022-30878-4), the corresponding authors Craig B. Wilen and Sidi Chen have listed their email addresses in their author information for that paper.

That being said, producing and testing an mRNA-LNP SARS-CoV-2 vaccine candidate is clearly within the capabilities of a small- to medium-scale team of microbiologists. It would be even easier if reviewers insisted that publications on the topic make a greater effort towards technology transfer, by doing things like publishing their construct sequences and the exact methodologies used to arrive at them.
