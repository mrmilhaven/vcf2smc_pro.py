#!/usr/bin/env python
import math
import numpy as np
import random
import vcf
import pysam
import sys
import time
names_input = argv[1]
vcf_file = 'NLE_noAsia.chr2.vcf.gz' #argv[2]
ancestral_bed = 'NLE_noAsia.chr2.bed' #argv[3]
converter = 'chr_converter.ancestralgibbon.txt'#argv[4]
ancestral_fa = '/vault/veeramah/gibbon/anc_seq/nomLeu3_ancestral.fa' #argv[5]
mask_file ='chr2.mask.bed'#argv[6]
pops = 'NLE_homozygous.RdistDip,NLE_homozygous.RdistHap'#argv[7]
output_stem = 'CUSTOM'#argv[8]

ref_ancestral = pysam.FastaFile(ancestral_fa)
vcf_reader = vcf.Reader(open(vcf_file,'r'))
temp_names = open(names_input,'r')
anc_bed = open(ancestral_bed,'r')
temp_converter = open(converter,'r')

begin_time = time.time()

individuals = []

class gibbon:
	population = 0
	name = ""
	ploidy = ""
	genotype = ""

u1_initial = 0
u2_initial = 0
global_alt_counter = 0
gib = gibbon()

pop=0
#Loads Names into List#
for line in temp_names:
	temp_pop = line.strip().split(" ")
	for item in temp_pop:
		name,ploidy = item.split(",")
		gib.name = name
		gib.ploidy = ploidy
		gib.population = pop

	
		if gib.population == 1:
			u1_initial = u1_initial+1 if gib.ploidy == "h" else u1_initial+2
			
		elif gib.population == 2:
			u2_initial = u2_initial+1 if gib.ploidy == "h" else u2_initial+2
		
		individuals.append(gib)
		gib = gibbon()
	pop += 1

chr_dict = {}
for line in temp_converter:

	vcf_chr_name,fasta_chr_name = line.split(" ")
	fasta_chr_name = fasta_chr_name.strip()
	chr_dict.update({vcf_chr_name:fasta_chr_name})


def most_frequent(gts):

	zero_count = 0
	one_count = 0
	two_count = 0
	three_count = 0

	for item in gts:

		zero_count += item.genotype.count('0')
		one_count += item.genotype.count('1')
		two_count += item.genotype.count('2')
		three_count += item.genotype.count('3')
	mfreq = max(zero_count,one_count,two_count,three_count)
	if mfreq == zero_count:
		num = 0
	elif mfreq == one_count:
		num = 1
	elif mfreq == two_count:
		num = 2
	elif mfreq == three_count:
		num = 3

	return num


def pop_creator(pop_list,pop):
        line=""
        add=""

	for item in pop_list:
		if item.population == pop:
			add = "[\""+item.name+"\", 0], [\""+item.name+"\", 1]"
        		add = add+", "
            		line = line + add
                #elif item != first:
                #    add = add+","
                
            		add = ""
	line = line[:-2]
	return line

	
pop1,pop2 = pops.split(",")
output_1 = open(str(pop1)+".smc" ,'w')
output_2 = open(str(pop2)+".smc",'w')
output_12 = open(str(pop1)+"_"+str(pop2)+".smc",'w')
output_21 = open(str(pop2)+"_"+str(pop1)+".smc",'w')


smcpp_beginning = "# SMC++ {\"version\": \"1.15.2\", \"pids\": "
dist_ind_name_pop1 = individuals[0].name
dist_ind_name_pop2 = individuals[1].name


output_1.write(smcpp_beginning+"[\""+pop1+"\"], \"undist\": [["+pop_creator(individuals,1)+"]], \"dist\": [[[\""+dist_ind_name_pop1+"\", 0], [\""+dist_ind_name_pop1+"\", 1]]]}")
output_1.write("\n")

output_2.write(smcpp_beginning+"[\""+pop2+"\"], \"undist\": [["+pop_creator(individuals,2)+"]], \"dist\": [[[\""+dist_ind_name_pop2+"\", 0], [\""+dist_ind_name_pop2+"\", 1]]]}")
output_2.write("\n")

output_12.write(smcpp_beginning+"[\""+pop1+"\", \""+pop2+"\"], \"undist\": [["+pop_creator(individuals,1)+"] ["+pop_creator(individuals,2)+"]], \"dist\": [[[\""+dist_ind_name_pop1+"\", 0], [\""+dist_ind_name_pop1+"\", 1]], [[\""+dist_ind_name_pop2+"\", 0], [\""+dist_ind_name_pop2+"\", 1]]]}")
output_12.write("\n")

output_21.write(smcpp_beginning+"[\""+pop2+"\", \""+pop1+"\"], \"undist\": [["+pop_creator(individuals,2)+"] ["+pop_creator(individuals,1)+"]], \"dist\": [[[\""+dist_ind_name_pop2+"\", 0], [\""+dist_ind_name_pop2+"\", 1]], [[\""+dist_ind_name_pop1+"\", 0], [\""+dist_ind_name_pop2+"\", 1]]]}")
output_21.write("\n")

def allele_counter(ancestral_num,individuals,u1_initial,u2_initial):
	d1 = 0
	d2 = 0
	u1 = u1_initial
	u2 = u2_initial
	n1 = u1_initial
	n2 = u2_initial
	isD1=True
	
	for item in individuals:
		
		if item.population == 0:
			if item.ploidy == 'h' and isD1:
                                if item.genotype == "0/1":
                                        if random.random() < 0.5:
                                                d1 = 1

                                elif item.genotype.count(".") == 0:
                                        if item.genotype.count(str(ancestral_num)) == 2:
                                                d1 = 0
					else:
						d1 = 1
				else:
					d1 = -1
                                isD1 = False 
                        elif item.ploidy == "h" and not isD1:
				if item.genotype == "0/1":
                                        if random.random() < 0.5:
                                                d2 = 1

                                elif item.genotype.count(".") == 0:
                                        if item.genotype.count(str(ancestral_num)) == 2:
                                                d2 = 0
					else:
						d2 = 1
				else:
					d2 = -1
			elif item.ploidy == "d" and isD1:
                                d1 = 2 - individuals[0].genotype.count(str(ancestral_num)) - item.genotype.count(".")
				isD1 = False
			elif item.ploidy == "d" and not isD1:
                                d2 = 2 - individuals[1].genotype.count(".")
		if item.population == 1:
			if item.ploidy == 'h':
				if item.genotype == "0/1":
					if random.random() < 0.5:
						u1 = u1-1
					
				elif item.genotype.count(".") == 0:
					if item.genotype.count(str(ancestral_num)) == 2:
						u1 = u1-1
				else:
					n1 = n1-1
					u1 = u1-1	
			else:
				u1 = u1 - item.genotype.count(str(ancestral_num)) - item.genotype.count(".")
				n1 = n1 - item.genotype.count(".")	

		elif item.population == 2:
			if item.ploidy == 'h':
                                if item.genotype == "0/1":
					if random.random() < 0.5:
						print("genotype decreased")
						u2 = u2-1

                                elif item.genotype.count(".") == 0:
					if item.genotype.count(str(ancestral_num)) == 2:
						u2 = u2-1
                                else:
                                        n2 = n2-1
                                        u2 = u2-1
                        else:
                                u2 = u2 - item.genotype.count(str(ancestral_num)) - item.genotype.count(".")
                                n2 = n2 - item.genotype.count(".")
	
	return d1,u1,n1,d2,u2,n2

def SNP_writer(ancestral_allele,individuals,record,u1_initial,u2_initial):

        global global_alt_counter
        ref = record.REF
        alts = record.ALT
        alt_counter = 1
    
	for item in individuals:
		item.genotype = record.genotype(item.name)['GT']


		if ref == ancestral_allele:
			line_to_write = allele_counter(0,individuals,u1_initial,u2_initial)
		elif ancestral_allele == 'N':
			line_to_write = allele_counter(most_frequent(individuals),individuals,u1_initial,u2_initial)
		else:
			for item in alts:
				if item == ancestral_allele:
					line_to_write = allele_counter(alt_counter,individuals,u1_initial,u2_initial)
					break
				else:
					line_to_write = allele_counter(most_frequent(individuals),individuals,u1_initial,u2_initial)
				alt_counter+=1
			global_alt_counter+=1
			alt_counter = 1	
	span = "1"
        d1 = str(line_to_write[0])
        u1 = str(line_to_write[1])
        n1 = str(line_to_write[2])
        d2 = str(line_to_write[3])
        u2 = str(line_to_write[4])
        n2 = str(line_to_write[5])

        pop1_line = " "+d1+" "+u1+" "+n1
        pop2_line = " "+d2+" "+u2+" "+n2
        smc1_line = span + pop1_line
        smc2_line = span + pop2_line
        smc12_line = span + pop1_line + pop2_line
        smc21_line = span + pop2_line + pop1_line

        print("ref:" + str(record.REF) + " alts:" + str(record.ALT) + " ancestral:" + str(ancestral_allele))
        print("SNP line:" + smc12_line)
        for item in individuals:
               print(item.name+":"+item.genotype)
	
	return smc1_line,smc2_line,smc12_line,smc21_line



def anc_writer(pos2,pos1,u1_initial,u2_initial):
	
	u1_initial = str(u1_initial)
	u2_initial = str(u2_initial)

	span = str(pos1 - pos2)
	du = " 0 0 "

	smc1_line = span+du+u1_initial
	smc2_line = span+du+u2_initial
	smc12_line = span+du+u1_initial+du+u2_initial
	smc21_line = span+du+u1_initial+du+u2_initial

	print("writing anc span...")
	return smc1_line,smc2_line,smc12_line,smc21_line

def mask_writer(pos2,pos1):
	
	span = str(pos1 - pos2 + 1)
	du = " -1 0 0"
	smc1_line = span+du
	smc2_line = span+du
	smc12_line = span+du+du
	smc21_line = span+du+du
	
	print("writing mask...")
	return smc1_line,smc2_line,smc12_line,smc21_line
	

record = next(vcf_reader)

prev_pos = 1
total_miss_maps = 0
total_span = 0
counter = 0

end_mask_prev = 0
try:
        temp_mask = open(mask_file,'r')
        #raw_input("Mask Found! Continue?...")
        mask_interval = next(temp_mask)
        mask_chromosome,begin_mask,end_mask = mask_interval.split("\t")
        end_mask = int(end_mask.strip())
        begin_mask = int(begin_mask)
        mask_lines = mask_writer(begin_mask,end_mask)

except:
        temp_mask = "none"
        raw_input("Mask not found. Is this correct?")

if temp_mask != "none":

    for line in anc_bed:
        counter +=1	
        chromosome_ancestor,start_pos_ancestor,end_pos_ancestor,original_coord_gibbon = line.split('\t')
        start_pos_ancestor = int(start_pos_ancestor)
        end_pos_ancestor = int(end_pos_ancestor)
        original_chr,current_pos = original_coord_gibbon.split("_")
        current_pos = current_pos.strip()
        current_pos = int(current_pos)
        while record.POS != current_pos:
            record = next(vcf_reader)
        
        ancestral_allele = ref_ancestral.fetch(chr_dict[chromosome_ancestor],start_pos_ancestor,end_pos_ancestor)
        print(str(start_pos_ancestor))	
        SNP_lines = SNP_writer(ancestral_allele,individuals,record,u1_initial,u2_initial)

        if current_pos > end_mask_prev and current_pos < begin_mask:
            
            anc_span = current_pos - prev_pos

            if anc_span > 0:
            
                anc_lines = anc_writer(prev_pos,current_pos,u1_initial,u2_initial)
                ###write anc###	
                output_1.write(anc_lines[0]+'\n')
                output_2.write(anc_lines[1]+'\n')
                output_12.write(anc_lines[2]+'\n')
                output_21.write(anc_lines[3]+'\n')
            ###WRITE SNP###	
            output_1.write(SNP_lines[0]+'\n')
            output_2.write(SNP_lines[1]+'\n')
            output_12.write(SNP_lines[2]+'\n')
            output_21.write(SNP_lines[3]+'\n')

            print("writing SNP...")
        elif current_pos > end_mask:

            anc_span = begin_mask - prev_pos
            
            if anc_span > 0:
            
                anc_lines = anc_writer(prev_pos,begin_mask,u1_initial,u2_initial)
                
                output_1.write(anc_lines[0]+'\n')
                output_2.write(anc_lines[1]+'\n')
                output_12.write(anc_lines[2]+'\n')
                output_21.write(anc_lines[3]+'\n')

            ###write Mask###
            output_1.write(mask_lines[0]+'\n')
            output_2.write(mask_lines[1]+'\n')
            output_12.write(mask_lines[2]+'\n')
            output_21.write(mask_lines[3]+'\n')

            ###NEW MASK###
            end_mask_prev = end_mask
            mask_interval = next(temp_mask)
            mask_chromosome,begin_mask,end_mask = mask_interval.split("\t")
            end_mask = int(end_mask.strip())
            begin_mask = int(begin_mask)
            mask_lines = mask_writer(begin_mask,end_mask)
            while current_pos > end_mask:
            
                anc_span = begin_mask - end_mask_prev
                
                if anc_span > 0:
                
                    anc_lines = anc_writer(end_mask_prev,begin_mask,u1_initial,u2_initial)
                    
                    output_1.write(anc_lines[0]+'\n')
                    output_2.write(anc_lines[1]+'\n')
                    output_12.write(anc_lines[2]+'\n')
                    output_21.write(anc_lines[3]+'\n')
            
                ###WRITE MASK###
                output_1.write(mask_lines[0]+'\n')
                output_2.write(mask_lines[1]+'\n')
                output_12.write(mask_lines[2]+'\n')
                output_21.write(mask_lines[3]+'\n')

                ###NEW MASK###
                end_mask_prev = end_mask
                mask_interval = next(temp_mask)
                mask_chromosome,begin_mask,end_mask = mask_interval.split("\t")
                end_mask = int(end_mask.strip())
                begin_mask = int(begin_mask)
                mask_lines = mask_writer(begin_mask,end_mask)
            anc_span = current_pos - end_mask_prev
                
            if anc_span > 0:
                anc_lines = anc_writer(end_mask_prev,current_pos,u1_initial,u2_initial)
                ###write anc###
                output_1.write(anc_lines[0]+'\n')
                output_2.write(anc_lines[1]+'\n')
                output_12.write(anc_lines[2]+'\n')
                output_21.write(anc_lines[3]+'\n')

            ###write SNP###
            output_1.write(SNP_lines[0]+'\n')
            output_2.write(SNP_lines[1]+'\n')
            output_12.write(SNP_lines[2]+'\n')
            output_21.write(SNP_lines[3]+'\n')
            
            print("writing SNP...")	
        prev_pos = current_pos
else:
    prev_pos = 0
    for line in anc_bed:
        counter +=1	
        chromosome_ancestor,start_pos_ancestor,end_pos_ancestor,original_coord_gibbon = line.split('\t')
        start_pos_ancestor = int(start_pos_ancestor)
        end_pos_ancestor = int(end_pos_ancestor)
        original_chr,current_pos = original_coord_gibbon.split("_")
        current_pos = current_pos.strip()
        current_pos = int(current_pos)
        print(current_pos,prev_pos)
        while record.POS != current_pos:
            record = next(vcf_reader)
        ancestral_allele = ref_ancestral.fetch(chr_dict[chromosome_ancestor],start_pos_ancestor,end_pos_ancestor)
        SNP_lines = SNP_writer(ancestral_allele,individuals,record,u1_initial,u2_initial)
        
        anc_span = current_pos - prev_pos
        
        if anc_span > 0:
            anc_lines = anc_writer(prev_pos,current_pos,u1_initial,u2_initial)
            ###write anc###	
            output_1.write(anc_lines[0]+'\n')
            output_2.write(anc_lines[1]+'\n')
            output_12.write(anc_lines[2]+'\n')
            output_21.write(anc_lines[3]+'\n')
        ###WRITE SNP###	
        output_1.write(SNP_lines[0]+'\n')
        output_2.write(SNP_lines[1]+'\n')
        output_12.write(SNP_lines[2]+'\n')
        output_21.write(SNP_lines[3]+'\n')
        prev_pos = current_pos
end_time = time.time()

print("the ancestral was not the reference "+str(global_alt_counter)+" times")
print("it took "+str(end_time - begin_time)+" seconds")	
