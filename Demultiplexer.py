# -*- coding: utf-8 -*-
"""
Created on April 2, 2019
@author: veej
File for demultiplexing and trimming the sequences.

In this script fastq files are analysed. This files contain out of DNA sequences, being 4 lines for 1 item.
> 1st line, consist out of a header.
> 2nd line, consist out of the DNA or RNA string-> this is the one we are interesed in and using in this script
> 3th line, only contains a '+' stupid useless and annoying I know..
> 4th line, has a quality score, also not used.
There is a R1 and a R2 file. They are complementary to each other. The first 4 lines of R1 are complementary
to the first 4 from R2.

###############################################################
# GENERAL WORKFLOW                                            #
#  0: Getting the files.                                      #
#  1: Print time and memory, to control duration + memory use #
#  2: Get the primers and barcodes from the mapfile.          #
#  3: create regex for the fw and rv primer + RC              #
#  4: Get 4 lines out of file, create a generator of them     #
#  5: Control sequences on N                                  #
#  6: Find primers in sequences                               #
#  7: Find barcode in sequences                               #
#  8: Trim barcode and primer of sequence                     #
#  9: Write the sequence in the correct file, by there status #
#  >: Repeat until all lines of the file are processed        #
# 10: Close al the used files                                 #
###############################################################

VERSION = 3.0 - 30 jun 2019
"""
# imports
# > Regex is used to search a small string inside a string
# > islice is used to slice the file in lines of 4
# > csv is used to open the map files as csv files
# > datetime is used to see the duration of the script
# > psutil is used to see the memory use of the script
# > glob is used to get the file names by pattern
# > os is used for removing the path names when writing the files and creating a directory
# > bioseq is used to get RC.
import regex as regex
from itertools import islice
import csv
import datetime
import psutil
import glob
import os
from Bio.Seq import Seq


# Do not modify

# constants to set log on or test on.
LOGGING_ENABLED = False
# setting normal and RC counter for the primer search
normal = 0
reverseComplement = 0
# setting the variables for the sequence status.
SEQUENCE_FORWARD = "forward"
SEQUENCE_FORWARD_RC = 'forward RC'
SEQUENCE_REVERSE_RC = 'reverse RC'
SEQUENCE_REVERSE = "reverse"
PRIMER_MISTAKE = "mistake"
BARCODE_MISTAKE = 'barcode'
N_DISCARD = "N"
# # set the number of mismatches
# num_mismatch = "1"
# print memory before
process = psutil.Process(os.getpid())
bef = process.memory_info().rss

def main():
    # print time and mem use
    print(datetime.datetime.now().time())
    print(bef)
    # get files
    R1_uri, R2_uri,num_mismatch, path, MAP_FILE,sample = collect_files()
    # get barcode and primers
    BARCODE_ID, fwPrimer, rvPrimer, fwrevcomprimer, rvrevcomprimer, file_dict, RC_BARCODE_ID = barcodes(MAP_FILE, path, sample)

    # make RE of the primers
    forward_primer_regex = makePrimerRE(fwPrimer,num_mismatch)
    reverse_primer_regex = makePrimerRE(rvPrimer,num_mismatch)
    fw_RC_primer_regex = makePrimerRE(fwrevcomprimer,num_mismatch)
    rv_RC_primer_regex = makePrimerRE(rvrevcomprimer,num_mismatch)

    # get the sequence
    sequences = getSequence(R1_uri, R2_uri)
    # filter the sequences with N out
    sequences = filter(lambda x: controlForN(x[0], x[1]), sequences)
    # search for primers
    primerMatches = (findPrimer(x[0], x[1], forward_primer_regex, reverse_primer_regex, fw_RC_primer_regex, rv_RC_primer_regex) for x in sequences)
    # search for barcodes
    barcodecheck = (findBarcode(x[1],x[2], BARCODE_ID, x[0], x[3], x[4], RC_BARCODE_ID) for x in primerMatches if x is not None)
    # trim the sequences
    trimSequences = (trimSeq(x[0], x[1], x[2], x[3], x[4], x[5]) for x in barcodecheck if x is not None)
    # control the status en write the sequence
    for i in trimSequences:
        controlSeqStatus(i[0],i[1],i[2], i[3], file_dict)
    # close all the used files
    closer(file_dict, path)


"""
In this function The R1, R2 FASTQ files and the mapfile is stored in a variable.
The sample names is created, being the first characters of the map file. 
Also the ammount of mismatches is set in a variable. Making it possible to specify this to a user input in a later version.
A directory is made for the output files if there is non.
three global variables are made for the mistake files.
"""
def collect_files():
    # setting the files, sample name and number of mismatches.
    R1_URI = glob.glob('*_R1_001.fastq')[0]
    R2_URI = glob.glob('*_R2_001.fastq')[0]
    MAP_FILE = glob.glob ('*.csv')[0]
    sample = MAP_FILE[:-12]
    num_mismatch = '1'

    # Creation of output dir. if there is none.
    if not os.path.isdir('./'+sample+'_output'):
        os.mkdir('./'+sample+'_output')
    path = './'+sample+'_output/'

    # For sequences with no primer, or only one primer find in the R1 but not in the complementary sequence
    global PRIMER_ERROR_FILE
    PRIMER_ERROR_FILE = open(path + "A3-Primer-Errors.fastq", "w")

    # For seqeunces with no barcode or two non overlapping barcodes.
    global BARCODE_ERROR_FILE
    BARCODE_ERROR_FILE = open(path+'A3-Barcode-Errors.fastq', 'w')

    # For all discarded sequences, being those containing a 'N'
    global DISCARD_FILE
    DISCARD_FILE = open(path+"A3-N-Discard.fastq", 'w')

    return(R1_URI, R2_URI, str(num_mismatch), path, MAP_FILE,sample)

"""
The barcode with sample IDs collected from the map file and stored in a dictionary
The output files, four per sample with all the possible sample and oriantation files 
are created and the paths are stored as well to make sure we find them back at the end.
The names are stored to control if there are actually samples inside at the end 
and to note which samples are non represented
The primers also collected out of the map file. 
The RC of the primers was made. 
"""
def barcodes(MAP_FILE,path,sample):
    with open(MAP_FILE, 'rU') as f:
        BARCODE_ID, RC_BARCODE_ID, file_dict= {}, {}, {}
        next(f)  # skip headers
        reader = csv.reader(f, delimiter='\t')
        # collectes barcode and sample ID, createas a dictionary with {barcode : sampleID}
        # also makes a RC dictionary of the barcodes.
        for sampleID, BarcodeSequence, LinkerPrimerSequence, ReversePrimer, *_ in reader:
            BARCODE_ID[BarcodeSequence] = sampleID
            RC_BARCODE_ID[str(Seq(BarcodeSequence).reverse_complement())] = sampleID
            # opens a four files per sample, one for every possible orientation
            r1fw = open(path + sampleID + '-R1fw.fastq', 'w+')
            r2rv = open(path + sampleID + '-R2rv.fastq', 'w+')
            r1rv = open(path + sampleID + '-R1rv.fastq', 'w+')
            r2fw = open(path + sampleID + '-R2fw.fastq', 'w+')
            # save the file names
            r1fname = sampleID + '-R1fw.fastq'
            r2rname = sampleID + '-R2rv.fastq'
            r1rname = sampleID + '-R1rv.fastq'
            r2fname = sampleID + '-R2fw.fastq'
            # make a list of the file names
            x = [r1fname, r2rname, r1rname, r2fname]
            # make a dictionary of the sample ID with all the files, and the list of file names.
            file_dict[sampleID] = [r1fw,r2rv,r1rv,r2fw,x]
        # For the SOIL FENCE data the fw and rv primer are selected from the mapfile from the 3th base,
        # Since the first 2 bases of the primer indicate RNA or DNA
        if 'SOIL_FENCE' in sample:
            fwPrimer = LinkerPrimerSequence[1:]
            rvPrimer = ReversePrimer[1:]
        # For the other datasets the primers are collected from the mapfile without modification
        else:
            fwPrimer = LinkerPrimerSequence
            rvPrimer = ReversePrimer
        # reverse complement is made of the primers
        fwrevcomprimer = Seq(fwPrimer).reverse_complement()
        rvrevcomprimer = Seq(rvPrimer).reverse_complement()

    return(BARCODE_ID, fwPrimer, rvPrimer, fwrevcomprimer, rvrevcomprimer, file_dict, RC_BARCODE_ID)

"""
From all the four possible primers a RE was created. 
Also enabling the possibility to process primers with e.g. R or Y in them, meaning that multiple nucleotide can 
be on one specific base. The RE was set to allow one mismatch. 
"""
def makePrimerRE(primer, num_mismatch):
    #start of RE
    rePrimer = "("
    # so I thought lets do a for loop, for every element out of the primer, if it is an A,T,G or C just make it a ATCG
    for i in primer:
        if i == "A":
            rePrimer += i
        elif i == "C":
            rePrimer += i
        elif i == "G":
            rePrimer += i
        elif i == "T":
            rePrimer += i
        # if it is something else lets check if it is one of the logical ones for sequences that then creates the RE.
        elif i == "R":
            rePrimer += "[AG]"
        elif i == "Y":
            rePrimer += "[CT]"
        elif i == "M":
            rePrimer += "[AC]"
        elif i == "K":
            rePrimer += "[GT]"
        elif i == "S":
            rePrimer += "[CG]"
        elif i == "W":
            rePrimer += "[AT]"
        elif i == "B":
            rePrimer += "[CGT]"
        elif i == "D":
            rePrimer += "[AGT]"
        elif i == "H":
            rePrimer += "[ACT]"
        elif i == "V":
            rePrimer += "[ACG]"
        # If it is something else it let us know.
        else:
            print("Unknown character find in primer sequence: " + i)
    # The whole RE is made, allowing the set numer of mismatches.
    rePrimer = regex.compile(rePrimer + ")" + "{e<="+num_mismatch+"}")

    return(rePrimer)

"""
The script reads one paired-end sequence at the time, containing the R1 read and related R2 read together with 
the header and quality score and analyses them all together.
So this function opens the R1 and R2 file, gets 4 lines out of the file in list format and creates a generator of them to continue.
stops when the end of the file is reached or when less then 4 lines are available.
"""
def getSequence(R1_uri, R2_uri):
    with open(R1_uri, 'rU') as R1_file, open(R2_uri, 'rU') as R2_file:
        while True:
            try:
                # get 4 lines out of the file as a list
                seq_R1 = list(islice(R1_file, 4))
                seq_R2 = list(islice(R2_file, 4))
                # control if there are 4
                if(len(seq_R1) < 4 or len(seq_R2) < 4):
                    # End of File
                    break
                # create a generator with the two  listed sequences.
                yield(seq_R1, seq_R2)

            except EOFError:
                #End Of File
                break

"""
The second step was to filter on "N" in the reads. If there was one or more inside one of the reads, 
this whole sequence received the status "N-Discard". If it passed the filter, the sequence moved on to the next step.
"""
def controlForN(seq_R1, seq_R2):
    if 'N' in seq_R1[1] or 'N' in seq_R2[1]:
        #LOG('Discard, an N has been found.', seq_R1[1])
        controlSeqStatus(N_DISCARD, seq_R1, seq_R2,0,0)
        return False
    else:
        return True

"""
Here the primer was searched in the first fifty base pairs of the reads. 
The location of the primers was saved and was later on used for trimming the sequences. 
First, the R1 read was controlled for the forward primer; if this primer was found the location was saved. 
Then the R2 read was controlled for the reverse primer, if this reverse primer was also found the location
was saved and the sequences received the status "Sequence forward". If only a forward primer was 
found in the R1 read the status became"Sequence mistake". This procedure was repeated for the 
reverse primer, the RC of the forward primer, and the RC of the reverse primer. 
In case both RC of the primers were found the status became sequence "forward-" or "reverse RC". 
To limit the processing time of the script the first fifty sequences were analysed on the original primers 
as on the RC of the primers, saving the number of original and RC primers found. 
After fifty sequences the script looked for the most often occurring orientation and continued searching 
for only that orientation of primers.
If no primers were found in R1, the sequence status becomes "Sequence mistake". 
In this case the locations were set to zero.
"""
# TODO: The control if 50 sequences are processed can be more general instead of checking every time..
# TODO: Search for the primers can be more general..
def findPrimer(seq_R1, seq_R2, forward_primer_regex, reverse_primer_regex, fw_RC_primer_regex, rv_RC_primer_regex):
    # RE search for fw, rv, and RC of fw and rv in first 50 nucleotides
    fw_found = regex.search(forward_primer_regex, seq_R1[1][:50])
    rv_found = regex.search(reverse_primer_regex, seq_R1[1][:50])
    fw_RC_found = regex.search(fw_RC_primer_regex, seq_R1[1][:50])
    rv_RC_found = regex.search(rv_RC_primer_regex, seq_R1[1][:50])

    #controls if 50 sequences have been processed.
    if (globals()['normal']+globals()['reverseComplement']) > 50:
        # controls which orientation is found the most
        if globals()['normal']>globals()['reverseComplement']:
            # wannebe bit control set to 1, to search only for normal orientation
            control = 1
        # if RC is found more, wannebe bit control is set to 2, to search only for the RC orientation
        else:
            control = 2
    # Keep going until 50 sequences are processed
    else:
        control = 0

    # controls for Fw in R1
    if fw_found != None and (control == 0 or control == 1):
        # save location
        locationR1 = fw_found.end()
        # search for rv in R2
        rv_also_found = regex.search(reverse_primer_regex, seq_R2[1][:50])
        if rv_also_found != None:
            # save location
            locationR2 = rv_also_found.end()
            # add 1 to global normal found orientation
            globals()['normal']+=1
            return (SEQUENCE_FORWARD, seq_R1, seq_R2,locationR1, locationR2)
        else:
            # Only the fw is found in R1, status becomes sequences mistake with location 0
            controlSeqStatus(PRIMER_MISTAKE, seq_R1, seq_R2, 0,0)

    # control for rv in R1
    elif rv_found != None and (control == 0 or control == 1):
        # save location
        locationR1 = rv_found.end()
        # search for fw in R2
        fw_also_found = regex.search(forward_primer_regex, seq_R2[1][:50])
        if fw_also_found != None:
            # save location
            locationR2 = fw_also_found.end()
            # add 1 to global normal found orienation
            globals()['normal']+=1
            return (SEQUENCE_REVERSE, seq_R1, seq_R2, locationR1, locationR2)
        else:
            # only the rv is found in R1, status becomes sequence mistake with location 0
            controlSeqStatus(PRIMER_MISTAKE, seq_R1, seq_R2,0,0)

    ###### ---- NOW FOR REVERSE COMPLEMENT ---- ###### -> could make this more general!!
    # search for RC fw in R1
    elif fw_RC_found != None and (control == 0 or control == 2):
        # save location
        locationR1 = fw_RC_found.end()
        # search for rc rv in R2
        rv_RC_also_found = regex.search(rv_RC_primer_regex, seq_R2[1][:50])
        if rv_RC_also_found != None:
            # save location
            locationR2 = rv_RC_also_found.end()
            # add 1 to global reverseComplement found orientation
            globals()['reverseComplement'] +=1
            return(SEQUENCE_FORWARD_RC, seq_R1, seq_R2, locationR1, locationR2)
        else:
            # only RC fw found in R1, status becomes sequence mistake
            controlSeqStatus(PRIMER_MISTAKE, seq_R1, seq_R2, 0, 0)

    # search for rc rv in R1
    elif rv_RC_found != None and (control == 0 or control == 2):
        # save location
        locationR1 = rv_RC_found.end()
        # search for rc fw in R2
        fw_RC_also_found = regex.search(fw_RC_primer_regex, seq_R2[1][:50])
        if fw_RC_also_found != None:
            # save location
            locationR2 = fw_RC_also_found.end()
            # add 1 to global reverse complement found orienation
            globals()['reverseComplement'] += 1
            return(SEQUENCE_REVERSE_RC, seq_R1, seq_R2, locationR1, locationR2)
        else:
            # only rc rv found in R1, status becomes sequence mistake
            controlSeqStatus(PRIMER_MISTAKE, seq_R1, seq_R2, 0, 0)

    # nothing found, status becomes sequence mistake
    else:
        controlSeqStatus(PRIMER_MISTAKE, seq_R1, seq_R2, 0, 0)

"""
The fourth step was the search for a barcode. The script searched for all the given barcodes one by one 
at the start of the reads, being the length of the barcode plus two. If the barcode was found, the barcode 
and sample ID were added to the header of the sequence. In the case that a second barcode was found the largest 
of the two was saved since barcodes can sometimes overlap. (I know, a bit annoying). If more than one barcode, that did not overlap, 
or no barcode was found the status of the sequence received the status "Barcode mistake".
"""
def findBarcode(seq_R1, seq_R2, BARCODE_ID, status, locationR1, locationR2, RC_BARCODE_ID):
    # set a barcode counter, because barcodes can be substrings of each oter......
    found_barcodes = []
    # get the barcodes and the correct dictionary for this sequence
    # if a sequence is RC, the RC dictionary is used. this holds the RC of the barcodes because everything is flipped..
    if status == SEQUENCE_REVERSE_RC or status == SEQUENCE_FORWARD_RC:
        keys = RC_BARCODE_ID.keys()
        dict = RC_BARCODE_ID
        print(status)
        print(seq_R2[0])
    # if the sequence is 'normal' the original orientation of the barcode is used.
    elif status == SEQUENCE_FORWARD or status == SEQUENCE_REVERSE:
        dict = BARCODE_ID
        keys = BARCODE_ID.keys()
    # some neurotic control else because I always want to keep track of everything.
    else:
        print("not sure what status we got to get the barcode. Let's check myself: ",status)

    # search for every barcode at the start of the R1 and R2 read.
    for i in keys:
        if i in seq_R1[1][:len(i) + 2] and i in seq_R2[1][:len(i) + 2]:
            # if found, add the barcode to the found barcode list
            found_barcodes.append(i)
            # control if ther is already an other barcode found
            if len(found_barcodes) == 2:
                # control if they are substrings if the new one is a substring of the old one, we use the old one.
                # bigger is better e.g.-> ATCGT or ATCGTCA then ATCGTCA used.
                if i in found_barcodes[0]:
                    # add the correct sample ID and barcode to the header,
                    seq_R1[0] = seq_R1[0][:-1] + '+' + dict[found_barcodes[0]] + '+' + found_barcodes[0] + '\n'
                    seq_R2[0] = seq_R2[0][:-1] + '+' + dict[found_barcodes[0]] + '+' + found_barcodes[0] + '\n'
                    # make the found barcodes again a list of the 1 correct found barcode
                    found_barcodes = [found_barcodes[0]]
                # if the old one is a supstring of the new one we use the new one
                elif found_barcodes[0] in i:
                    # add the sample ID and barcode to the header
                    seq_R1[0] = seq_R1[0][:-1] + '+' + dict[i] + '+' + i + '\n'
                    seq_R2[0] = seq_R2[0][:-1] + '+' + dict[i] + '+' + i + '\n'
                    # set the found barcode list to the found barcode.
                    found_barcodes = [i]

                else:
                    # two barcode founds that do not overlap. The status becomes barcode mistake
                    controlSeqStatus(BARCODE_MISTAKE, seq_R1, seq_R2,0,0)
            # the first barcode found, so no difficult things. just add the sample ID and barcode to the header
            else:
                seq_R1[0] = seq_R1[0][:-1] + '+' + dict[i] + '+' + i + '\n'
                seq_R2[0] = seq_R2[0][:-1] + '+' + dict[i] + '+' + i + '\n'

    # control at the end if none or multiple not overlapping barcodes are found (last is impossible to have here..)
    # if so status is set to barcode mistace
    if len(found_barcodes) == 0 or len(found_barcodes) >= 2:
        controlSeqStatus(BARCODE_MISTAKE, seq_R1, seq_R2,0,0)
    else:
        return(status, seq_R1, seq_R2, locationR1, locationR2, dict[found_barcodes[0]])

"""
The fifth step in the program was the trimming of the sequences. The saved locations of the found primers 
were used to trim the sequences as well as the corresponding quality scores.
In the case of a RC primer being, the sequence was trimmed and the the RC was taken of the whole sequence. 
In consequence the quality score was reversed to ensure the correspondence of the quality score per base 
with the correct base in the read.
"""
# TODO: the whole RC part feels a bit odd...
def trimSeq(status, seq_R1_Barcode,seq_R2_Barcode,locationR1, locationR2, sampleID):
    # trim the sequences
    seq_R1_Barcode[1] = seq_R1_Barcode[1][locationR1:]
    seq_R2_Barcode[1] = seq_R2_Barcode[1][locationR2:]
    seq_R1_Barcode[3] = seq_R1_Barcode[3][locationR1:]
    seq_R2_Barcode[3] = seq_R2_Barcode[3][locationR2:]

    # in case of a RC, flip the sequences and quality scores.
    if status == SEQUENCE_FORWARD_RC or status == SEQUENCE_REVERSE_RC:
        seq_R1_Barcode[1] = Seq(seq_R1_Barcode[1]).reverse_complement()
        seq_R2_Barcode[1] = Seq(seq_R2_Barcode[1]).reverse_complement()
        seq_R1_Barcode[3] = seq_R1_Barcode[3][::-1]
        seq_R2_Barcode[3] = seq_R2_Barcode[3][::-1]

    return(status, seq_R1_Barcode,seq_R2_Barcode, sampleID)

"""
In the last step, the status of the sequences was controlled. From there on the right function was called
Because there were some problems with lengths of sequences and quality scores, there is a control for them here.
The sequences are only written if the read and quality score are of the same length.
Write the sequences in the appropriate orientation file. 
The sequences with a "N", primer mistake or barcode mistake were written in their own files. 
These mistake files have not been used in the ESV analysis but are used and stored to provide 
insight into the quality of the data set.
"""
def controlSeqStatus(seq_status, seq_R1, seq_R2, sampleID, file_dict):
    # control length read and qualtiy score
    if len(seq_R1[1].encode('ascii')) == len(seq_R1[3].encode('ascii')) and len(seq_R2[1].encode('ascii')) == len(seq_R2[3].encode('ascii')):
        # if length is the same, control the status and write the sequence in the correct file.
        # The forward seauences, writes the 4 sequences lines in the R1 forward and R2 reverse files.
        # This means that R1 has the forward primer.
        if seq_status == SEQUENCE_FORWARD:
            r1fwfile = file_dict[sampleID][0]
            r2rvfile = file_dict[sampleID][1]
            r1fwfile.writelines(seq_R1)
            r2rvfile.writelines(seq_R2)

        # The reverse sequences, writes the 4 sequences lines in the R1 revers and R2 forward files.
        # This means that R1 has the reverse primer.
        elif seq_status == SEQUENCE_REVERSE:
            r1rvfile = file_dict[sampleID][2]
            r2fwfile = file_dict[sampleID][3]
            r1rvfile.writelines(seq_R1)
            r2fwfile.writelines(seq_R2)

        # The primer mistake sequences. This one is called when no primer is found in the sequence.
        # It writes the 4 sequences lines in the ERROR file.
        elif seq_status == PRIMER_MISTAKE:
            PRIMER_ERROR_FILE.writelines(seq_R1)
            PRIMER_ERROR_FILE.writelines(seq_R2)

        # The barcode mistake sequences. This one is called when two or no barcode is found in the sequence.
        # It writes the 4 sequences lines in the barcode error file.
        elif seq_status == BARCODE_MISTAKE:
            BARCODE_ERROR_FILE.writelines(seq_R1)
            BARCODE_ERROR_FILE.writelines(seq_R2)

        # The discard sequences. This one is called when a N is found in the sequence or an other problem is found.
        # It writes the 4 sequences lines in the discard file.
        elif seq_status == N_DISCARD:
            DISCARD_FILE.writelines(seq_R1)
            DISCARD_FILE.writelines(seq_R2)
    # if the length of the read and quality is not the same, the sequence is written in the discard file.
    else:
        DISCARD_FILE.writelines(seq_R1)
        DISCARD_FILE.writelines(seq_R2)

"""
Close all the files :)
it also deletes all the empty files, since all possibilities are made and not all possibilities are used.
"""
def closer(file_dict, path):
    PRIMER_ERROR_FILE.close()
    BARCODE_ERROR_FILE.close()
    DISCARD_FILE.close()
    # make a file to write down the not represented samples
    with open(path + 'A3-Empty.fastq', 'w+') as e:
        # control if there is a seq. in every file, if so close it.
        # if not delete it and write the file name in the error file
        for r1fw,r2rv,r1rv,r2fw,names in file_dict.values():
            if r1fw.truncate() == 0:
                e.write(names[0] + '\n')
                os.remove(path +names[0])
                r1fw.close()
            if r2rv.truncate() == 0:
                e.write(names[1] + '\n')
                os.remove(path + names[1])
                r2rv.close()
            if r1rv.truncate() == 0:
                e.write(names[2] + '\n')
                os.remove(path + names[2])
                r1rv.close()
            if r2fw.truncate() == 0:
                e.write(names[3] + '\n')
                os.remove(path + names[3])
                r2fw.close()
            else:
                r1fw.close()
                r2rv.close()
                r1rv.close()
                r2fw.close()

    print(datetime.datetime.now().time())
    af = process.memory_info().rss
    print(af)
    print(af-bef)

# function to print logs for on the flow control of anything.
# e.g.  put LOG("this will be printed if log enabled") anywhere to log something during programming.
def LOG(*messages):
    if LOGGING_ENABLED:
        print(messages)

# The always forgotten call of the main function.
main()
