#!/usr/bin/gawk -f

# Script to remove reads from a sam file with indels
# Dr. Sreenu Vattipally
# Bioinformatics Hub
# MRC-CVR, University of Glasgow Centre for Virus Research

#
# Usage: SamRemoveIndels.awk file.sam
#

# Field separator is a TAB
BEGIN { FS="\t"; }

# Print header lines
$1~/^@/{ print $0}

# Not a header and Mapped
$1!~/^@/ && !and($2,0x4){

	# If forward read
	if(!and($2,0x10)){

		read[$1,1]=$0;
		cigar[$1,1]=$6;
	}
	# If reverse read
	else if(and($2,0x10)){

		read[$1,2]=$0;
		cigar[$1,2]=$6;
	}
}

END{
	for(j in read){ 
		split(j,name,SUBSEP)
		if(cigar[name[1],1]!~/[ID]/ && cigar[name[1],2]!~/[ID]/){
			if(read[name[1],1]!="") {printf("%s\n",read[name[1],1]); read[name[1],1]="";}
			if(read[name[1],2]!="") {printf("%s\n",read[name[1],2]); read[name[1],2]="";}
		}
	}
}
