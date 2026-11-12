## Making the lepr exons bed files ##

Check recent R149 generate bed job, the transcript used for LEPR is NM_002303.6

**GRCh37_LEPR_exons_plus25.bed**
- Get the regions from GCF_000001405.25_GRCh37.p13_genomic_20201022.symbols.exon_5bp.tsv
- dx cat file-J1gqjQQ49jjBgzJQYzvG89YB | grep NM_002303.6 > GRCh37_LEPR_exons_plus25.bed
- -20 bases from start coordinate, +20 bases to end coordinate

**GRCh38_LEPR_exons_plus25.bed**
- Get the regions from GCF_000001405.39_GRCh38.p13_genomic_20211119.exon_5bp.tsv
- dx cat file-GyFfgpQ4fJPv132574bFQfV5 | grep NM_002303.6 > GRCh38_LEPR_exons_plus25.bed
- -20 bases from start coordinate, +20 bases to end coordinate