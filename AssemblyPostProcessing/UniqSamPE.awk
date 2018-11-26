#!/usr/bin/gawk -f

# Script to remove duplicate entries from a sam file
# Dr. Sreenu Vattipally
# Bioinformatics Hub
# MRC-University of Glasgow Centre for Virus Research

#
# Usage: UniqSamPE.awk file.sam
#

# Field separator is a TAB
BEGIN { FS="\t"; }

# Print header lines
$1~/^@/{ print $0}

# Not a header and Mapped
$1!~/^@/ && !and($2,0x4){

	# If forward read
	if(!and($2,0x10)){
		# If first in the read
		#if(and($2,0x40)) names[$1,1]=$4;
		# If second in the read
		#if(and($2,0x40)) names[$1,1]=$4 + length($10);

		names[$1,1]=$4;
		entry[$1,1]=$0;
	}
	# If reverse read
	else if(and($2,0x10)){
		# If first in the read
		#if(and($2,0x40)) names[$1,2]=$4;
		# If second in the read
		#else if(and($2,0x80)) names[$1,2]=$4+length($10);

		names[$1,2]=$4+length($10);
		entry[$1,2]=$0;
	}
}

END{
	for(i in names){
		split(i,indx,SUBSEP)
		position[names[indx[1],1],names[indx[1],2]]=1;
	}
	for(j in names){ 
		split(j,name,SUBSEP)
		if(position[names[name[1],1],names[name[1],2]]==1){
			if(entry[name[1],1]!="")
			printf("%s\n",entry[name[1],1]);
			if(entry[name[1],2]!="")
			printf("%s\n",entry[name[1],2]);
			position[names[name[1],1],names[name[1],2]]=0;
			position[names[name[1],2],names[name[1],1]]=0;
		}
	}
}
