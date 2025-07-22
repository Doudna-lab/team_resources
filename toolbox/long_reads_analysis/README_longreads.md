# Long Reads Analysis

## Data Acquisition

+ Oxford Nanopore Data was acquired via lftp
```angular2html
lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u o250714_Guerrero,iRie9eid6goHaik -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'
```