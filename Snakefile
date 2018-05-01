#!/usr/bin/env python3

import os


#############
# FUNCTIONS #
#############

def find_fastq_files(fastq_dir):
    fastq_dir_files = list((dirpath, filenames)
                           for (dirpath, dirnames, filenames)
                           in os.walk(fastq_dir))
    my_fastq_files = []
    for dirpath, filenames in fastq_dir_files:
        for filename in filenames:
            if filename.endswith('.fastq'):
                my_fastq_files.append(os.path.join(dirpath, filename))
    return(sorted(set(my_fastq_files)))


###########
# GLOBALS #
###########

sample_list = ['sample_3',
               'sample_5',
               'sample_7',
               'sample_8',
               'sample_10',
               'blue_cod']

compare_list = ['lambda_qc',
                'mh_78',
                'asw_11',
                'asw_12b',
                'asw_12a',
                'asw_14',
                'asw_31',
                'asw_33',
                'asw_34',
                'bm_c',
                'bm_wt',
                'lug_8']

#########
# RULES #
#########

rule target:
    input:
        'output/020_stats/combined/hist.pdf',
        'output/010_minion-qc/combined/qc_stats.png',
        expand('output/020_stats/{sample}/{type}_hist.pdf',
               sample=sorted(set(sample_list + compare_list)),
               type=['all', 'pass'])

rule plot_combined_hist:
    input:
        lhist = expand('output/020_stats/{sample}/{type}_hist.txt',
                       sample=sorted(set(sample_list + compare_list)),
                       type=['all', 'pass']),
    output:
        lh = 'output/020_stats/combined/hist.pdf',
        wh = 'output/020_stats/combined/weighted-hist.pdf',
        png = 'output/020_stats/combined/weighted-hist.png'
    log:
        log = 'output/logs/020_stats/combined_histogram.log'
    threads:
        1
    script:
        'src/plot_combined_histograms.R'

# plot histos, pass hist files and output files
rule plot_hist:
    input:
        lhist = 'output/020_stats/{sample}/{type}_hist.txt'
    output:
        lh = 'output/020_stats/{sample}/{type}_hist.pdf',
        wh = 'output/020_stats/{sample}/{type}_weighted-hist.pdf'
    log:
        log = 'output/logs/020_stats/{sample}-{type}_histogram.log'
    threads:
        1
    script:
        'src/plot_histograms.R'


rule all_stats:
    input:
        fail_dir = 'data/workshop2018/{sample}/basecalled/workspace/fail',
        pass_dir = 'data/workshop2018/{sample}/basecalled/workspace/pass'
    output:
        all_hist = 'output/020_stats/{sample}/all_hist.txt'
    threads:
        1
    log:
        'output/logs/020_stats/{sample}-pass_reformat.log'
    threads:
        1
    run:
        pass_files = find_fastq_files(input.pass_dir)
        fail_files = find_fastq_files(input.fail_dir)
        shell('cat {pass_files} {fail_files} | '
              'reformat.sh '
              'in=stdin.fastq '
              'interleaved=f '
              'out=stdout.fastq '
              'lhist={output.all_hist} '
              'maxhistlen=200000 '
              '> /dev/null '
              '2> {log}')


rule pass_stats:
    input:
        pass_dir = 'data/workshop2018/{sample}/basecalled/workspace/pass'
    output:
        pass_hist = 'output/020_stats/{sample}/pass_hist.txt'
    log:
        'output/logs/020_stats/{sample}-pass_reformat.log'
    threads:
        1
    run:
        pass_files = find_fastq_files(input.pass_dir)
        shell('cat {pass_files} | '
              'reformat.sh '
              'in=stdin.fastq '
              'interleaved=f '
              'out=stdout.fastq '
              'lhist={output.pass_hist} '
              'maxhistlen=200000 '
              '> /dev/null '
              '2> {log}')

rule parse_qc:
    input:
        yaml = expand('output/010_minion-qc/{sample}/summary.yaml',
                      sample=sorted(set(sample_list + compare_list)))
    output:
        parsed_data = 'output/010_minion-qc/combined/qc_data.Rds',
        gp = 'output/010_minion-qc/combined/qc_stats.pdf',
        png = 'output/010_minion-qc/combined/qc_stats.png'
    log:
        log = 'output/logs/010_minion-qc/parse_qc.log'
    threads:
        1
    script:
        'src/parse_minion_qc.R'

rule qc:
    input:
        seq_sum = 'data/workshop2018/{sample}/basecalled/sequencing_summary.txt'
    output:
        'output/010_minion-qc/{sample}/summary.yaml'
    params:
        outdir = 'output/010_minion-qc/{sample}'
    log:
        'output/logs/010_minion-qc/{sample}.log'
    shell:
        'Rscript src/MinIONQC.R '
        '-i {input.seq_sum} '
        '-o {params.outdir} '
        '&> {log}'
