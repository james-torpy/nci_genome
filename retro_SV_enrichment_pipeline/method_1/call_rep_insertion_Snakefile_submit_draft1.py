#!/usr/bin/env/python

import os
import re
import glob

project_dir = '/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/'
bwa_dir = project_dir + 'results/bwa/'
tape_dir = 'jt3341/projects/hgsoc_repeats/genome/raw_files/bowtell/'
svaba_dir = project_dir + 'results/svaba/'
cnvnator_dir = project_dir + 'results/cnvnator/'

# determine number of bams in bwa directory:
init_bam_no = len(glob.glob1(bwa_dir,'*.bam'))

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if not os.path.isfile(project_dir + 'done_list.txt'):
    os.system('echo done >>' + project_dir + 'done_list.txt')

# if 0 bams in directory, load next 15 jobs:
if init_bam_no == 0:
    print('\nNo bams left in in_dir - loading next 15 jobs...\n')
    global SAMPLES
    SAMPLES = []
    if file_len(project_dir + 'done_list.txt') < file_len(project_dir + 
        'bam_list.txt') + 1:
        with open(project_dir + 'bam_list.txt', 'r') as file_content:
            line = re.sub("\n", "", file_content.readline())
            while line:
                print('Sample name is: ' + str(line))
                if str(line) in open(project_dir + 'done_list.txt').read():
                    print(str(line) + ' has already been processed, moving on ' + 
                        'to next sample...')
                    line = re.sub("\n", "", file_content.readline()) 
                else:
                    # determine number of bams in bwa directory:
                    bam_no = len(glob.glob1(bwa_dir,'*.bam'))
                    if bam_no < 15:
                        print('Number of bams in in_dir: ' + str(bam_no) + 
                            '\n')
                        if not (os.path.isfile(svaba_dir + str(line) + 
                            '.discordant.txt.gz') and os.path.isfile(cnvnator_dir + str(line) + 
                            '.out.root')):
                            if not os.path.isfile(bwa_dir + str(line) + '.bam'):
                                # get bam file from tape:
                                bamfile = tape_dir + str(line)    + '.bam"'
                                print('Fetching ' + bamfile + ' from tape\n')
                                os.system('cd ' + bwa_dir + '; mdss get ' + tape_dir + 
                                    str(line)    + '.bam')
                            # append sample name to 'SAMPLES' array:
                            print('Appending ' + str(line)    + ' to SAMPLES array...\n\n')
                            SAMPLES.append(str(line)   )
                        else:
                            print(svaba_dir + str(line) + '.discordant.txt.gz' + 
                                ' and ' + cnvnator_dir + str(line) + '.out.root' + 
                                ' already exist\n\n')
                    else:
                        print('15 files already in directory\n\n')
    #                    files = glob.glob1(bwa_dir,'*.bam')
    #                    done_files = str(open(project_dir + 
    #                        'done_list.txt').readlines())
    #                    done_files = list(map(lambda x: x.replace('\n', ''), done_files))
    #                    for f in files:
    #                        sname = re.sub('.bam', '', f)
    #                        if not (os.path.isfile(svaba_dir + str(sname) + 
    #                            '.discordant.txt.gz') and os.path.isfile(cnvnator_dir + 
    #                            str(sname) + '.out.root') and str(sname) not in files):
    #                            SAMPLES.append(str(sname))
                    line = re.sub("\n", "", file_content.readline()) 
        # convert SAMPLES from list to string:
        SAMPLES = '_'.join(SAMPLES)
        # pass sample names to file:
        print('samples are ' + SAMPLES)
    
        SAMPLES_file = open(project_dir + 'SAMPLES.txt', 'w')
        SAMPLES_file.write(SAMPLES)
        SAMPLES_file.close()
        
        # submit snakefile for rep insert pipeline:
        os.system('cd ' + project_dir + '; qsub ' + project_dir + 
            'rep_insertion_Snakefile_submit.bash')
        
    