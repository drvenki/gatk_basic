#! /bin/bash
# User Input Look up file (vcf) 
read -p "Enter VCF file name:" lookup
# Gene of Interest - User Input
read -p "Type your gene of interest:" IN

while IFS=',' read -ra ADDR; do
      for i in "${ADDR[@]}"; do
      	echo -e "CHROM REF ALT GENE G.Syntax P.Syntax TYPE IMPACT"
        grep -v "^#" $lookup | \
        	awk -v gene=$i -F "\t" 'n=split($8,a,","){
        		for(i=1; i<=n; i++){ 
        			split(a[i],data,"|"); 
        			if (data[4]==gene) { 
        				gsub("ANN=","",data[1]); 
        				print $1,$4,$5,data[4],data[10],data[11],data[2],data[3] 
        			} 
        		} 
        	}'
      done 
done <<< "$IN"

