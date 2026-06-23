#! /usr/bin/env python3.7
# Copyright (C) 2024 Alvaro Chiner Oms & Miguel Moreno Molina

# This file is part of ThePipeline2
# ThePipeline2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# ThePipeline2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with ThePipeline2.  If not, see <http://www.gnu.org/licenses/>.

# Module for variant calling, using Varscan2, Mutect2 and Minos

import vcf
def check_extension(ext):
    '''Check if the extension is either BAM or CRAM'''
    import sys

    if ext.endswith(".bam") or ext.endswith(".cram"):
        return 0
    else:
        sys.exit("\033[91mERROR: We only accept files"
                 " ended in .bam or .cram extension\033[0m")


def VCFtoPandas(file):
    import pandas
    import io

    with open(file, 'r') as f:
        lines = [li for li in f if not li.startswith('##')]
    data = pandas.read_csv(io.StringIO(''.join(lines)),
                           dtype={'#CHROM': str, 'POS': int, 'ID': str,
                                  'REF': str, 'ALT': str, 'QUAL': str,
                                  'FILTER': str, 'INFO': str},
                           sep='\t')

    return data

def VarScan(reference, prefix, varscan, samtools, ext):
    '''Call samtools mpileup and Varscan, then transform output to simple VCF for Minos'''
    from subprocess import call
    import sys
    import pandas as pd
    from .History import UpdateHistory

    if ext == ".sort.bam":
        bamfile = "{}{}".format(prefix, ext)
        with open("{}.mpileup.remove".format(prefix), "w") as outfh:
            cmd = [samtools, "mpileup", "-AB", "-f", reference, bamfile]
            stat = call(cmd, stdout=outfh)
            if stat != 0:
                sys.exit("Pipeline stopped at samtools mpileup!\n")
            UpdateHistory(cmd, "samtools", prefix)

    # Si el archivo es .cram, lo convertimos a .sort.bam
    elif ext == ".cram":
        print("Converting {}.cram to {}.sort.bam".format(prefix,prefix))
        cmd = [samtools, "view", "-T", reference, "-B", "-o", "{}.sort.bam".format(prefix), "{}.cram".format(prefix)]
        call(cmd)
        bamfile = "{}{}".format(prefix, ext)
        with open("{}.mpileup.remove".format(prefix), "w") as outfh:
            cmd = [samtools, "mpileup", "-AB", "-f", reference, bamfile]
            stat = call(cmd, stdout=outfh)
            if stat != 0:
                sys.exit("Pipeline stopped at samtools mpileup!\n")
            UpdateHistory(cmd, "samtools", prefix)


    with open("{}.snp".format(prefix), "w") as outfh:
        cmd = ["java", "-Xms10G", "-Xmx32G", "-jar", varscan,
            "pileup2snp", "{}.mpileup.remove".format(prefix),
            "--min-coverage", "3",
             "--min-reads2", "3",
             "--min-freq-for-hom", "0.90",
             "--min-var-freq", "0.05"]
        stat = call(cmd, stdout=outfh)
        if stat != 0:
            sys.exit("Pipeline stopped at VarScan pileup2snp!\n")
        UpdateHistory(cmd, "varscan", prefix)
        call(["rm", "{}.mpileup.remove".format(prefix)])

        UpdateHistory(cmd, "varscan", prefix)

    with open("{}.parsed.vcf".format(prefix), "w") as outfh:
        outfh.write("##fileformat=VCFv4.2\n")
        outfh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        outfh.write("##contig=<ID=MTB_anc,length=4411532,assembly=MTB_ancestor_reference.fasta>\n")
        outfh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + prefix + "\n")
        snpfile = pd.read_csv("{}.snp".format(prefix), sep='\t', header=0)
        for index, snp in snpfile.iterrows():
            depth = str(int(snp["Reads1"]) + int(snp["Reads2"]))
            f = str(float(snp["VarFreq"].replace("%","").replace(",","."))/100)
            outfh.write("\t".join([snp["Chrom"], str(snp["Position"]), ".", snp["Ref"],
                        snp["VarAllele"], ".", ".", ".", "GT:AD:AF", "0/1:" + depth + ":" + f + "\n"]))


def Mutect2(reference, prefix, gatk, samtools, genomeCoverageBed,
            ext, threads, min_depth, min_qual, freq):
    '''Call Mutect2 to perfom variant calling and obtain a VCF file
    and bedtools to calculate coverage and determine lowcov positions'''
    import subprocess as sp
    import sys
    import os
    import pandas
    from .History import UpdateHistory

    # need to import libstdc++.so.6 updated library for
    # bedtools to work. It's in /PipeModules/Configs/lib
    dirname = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
    os.environ['LD_LIBRARY_PATH'] = os.path.join(dirname,
                                                 'data',
                                                 'libs')

    rm_bam = False
    if ext.endswith(".cram"):  # if input is CRAM, transform to BAM
        cmd_toBam = [samtools, "view", "--threads", threads,
                     "-T", reference,
                     "-B", "-o", "{}.sort.bam".format(prefix),
                     "{}{}".format(prefix, ext)]
        stat = sp.call(cmd_toBam)
        UpdateHistory(cmd_toBam, "samtools", prefix)
        if stat != 0:
            sys.exit("\033[91mPIPELINE ERROR at sample {}: "
                     "Pipeline stopped at samtools "
                     "when transforming BAM to CRAM!\n\033[0m".format(prefix))
        ext = ".sort.bam"

    # index BAM
    cmd_indexBam = [samtools, "index", "-@", threads,
                    "{}{}".format(prefix, ext)]
    stat = sp.call(cmd_indexBam)
    UpdateHistory(cmd_indexBam, "samtools", prefix)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: "
                 "Pipeline stopped when indexing BAM!\n\033[0m".format(prefix))

    # CALLING WITH MUTECT2
    # Run Mutect2
    cmd_Mutect = [gatk, "Mutect2",
                   "-R", reference, "-I", "{}{}".format(prefix, ext),
                   "-O", "{}.vcf".format(prefix), "-OVI", "false",
                   "--verbosity", "ERROR",
                   "--QUIET", "true", "-mbq", min_qual,
                   "--callable-depth", min_depth,
                   "-mnp-dist", "0",  # important to not filter MNPs later
                   "--linked-de-bruijn-graph", "true",
                   "--native-pair-hmm-threads", threads]
    
    stat = sp.call(cmd_Mutect)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: Pipeline stopped"
                 " when running Mutect2!\n\033[0m".format(prefix))
    UpdateHistory(cmd_Mutect, "gatk", prefix)
    
    #now produce the gVCF for accurate WT calling
    print("\033[92m\nObtaining {}.gvcf\n\033[00m".format(prefix))
    cmd_Mutect2 = [gatk, "Mutect2",
                   "-R", reference, "-I", "{}{}".format(prefix, ext),
                   "-O", "{}.gvcf".format(prefix), "-OVI", "false",
                   "--verbosity", "ERROR",
                   "--QUIET", "true", "-mbq", min_qual,
                   "--callable-depth", min_depth,
                   "-mnp-dist", "0",  # important to not filter MNPs later
                   "--linked-de-bruijn-graph", "true",
                   "--native-pair-hmm-threads", threads,
                   "-ERC", "BP_RESOLUTION"]
    
    stat = sp.call(cmd_Mutect2)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: Pipeline stopped"
                 " when running Mutect2 for gVCF generation!\n\033[0m".format(prefix))
    UpdateHistory(cmd_Mutect2, "gatk", prefix)

    # Run SelectVariants to split INDELs and SNPs
    cmd_selectI = [gatk, "SelectVariants", "-R", reference,
                   "-V", "{}.vcf".format(prefix),
                   "--verbosity", "ERROR",
                   "--QUIET", "true",
                   "-OVI", "false",
                   "-select-type", "INDEL", "-O",
                   "{}.indel.vcf".format(prefix)]
    stat = sp.call(cmd_selectI)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: "
                 "Pipeline stopped at SelectVariants\n\033[0m".format(prefix))
    UpdateHistory(cmd_selectI, "gatk", prefix)

    cmd_selectS = [gatk, "SelectVariants", "-R", reference,
                   "-V", "{}.vcf".format(prefix),
                   "--verbosity", "ERROR",
                   "--QUIET", "true",
                   "-OVI", "false",
                   "-xl-select-type", "INDEL",
                   "-xl-select-type", "MNP",
                   "-xl-select-type", "SYMBOLIC",
                   "-xl-select-type", "NO_VARIATION",
                   "-O", "{}.snp.vcf".format(prefix)]
    stat = sp.call(cmd_selectS)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: "
                 "Pipeline stopped at SelectVariants\n\033[0m".format(prefix))
    UpdateHistory(cmd_selectS, "gatk", prefix)

    # Extract multiallelic from indel.vcf
    sp.run("grep -v 0/1/2 {0}.indel.vcf > {0}.indel_sin_multiallelic.vcf".format(prefix), shell=True)
    UpdateHistory("grep -v 0/1/2 {0}.indel.vcf > {0}.indel_sin_multiallelic.vcf".format(prefix), "custom", prefix)
    # Filter indels 
    indel_reader = vcf.Reader(open("{}.indel_sin_multiallelic.vcf".format(prefix), "r"))
    indel_vcf_remade = vcf.Writer(open("{}.remade.indel.vcf".format(prefix), "w+"), indel_reader)

    for record in indel_reader:
        if(record.samples[0]['AF'] >= 0.05 and record.samples[0]['DP'] >= 10):
            indel_vcf_remade.write_record(record)

    indel_vcf_remade.close()

    # CALCULATE LOWCOV FROM COVERAGE
    # Calculate coverage
    with open("{}.coverage".format(prefix), "w") as outfh:
        cmd_cov = [genomeCoverageBed, "-ibam", "{}{}".format(prefix, ext),
                   "-d", "-g", reference]
        stat = sp.call(cmd_cov, stdout=outfh)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: Pipeline stopped when"
                 " calculating BAM coverage!\n\033[0m".format(prefix))
    outfh.close()
    UpdateHistory(cmd_cov, "bedtools", prefix)

    # lowcov will be all bases bellow min_depth cut-off
    cov_info = pandas.read_csv("{}.coverage".format(prefix),
                               sep='\t')
    cov_info.columns = ['CHROM', 'POS', 'COV']
    cov_info[cov_info['COV'] < float(min_depth)
             ]['POS'].to_csv("{}.lowcov.tsv".format(prefix),
                             sep='\n',
                             index=False)


# Filter multiallelic from Mutect2 snp output
def filter_multiallelic_from_mutect2_snp(prefix):
    from .History import UpdateHistory
    '''For every multiallelic position in .snp.vcf from Mutect2 keeps the SNP with the highest freq value'''

    with open("{}.snp.vcf".format(prefix), "r+") as input_file:
        lines = input_file.readlines()
    
    headers = [line for line in lines if line.startswith("#")]
    lines_header_CHROM = [index for index, element in enumerate(headers) if "#CHROM" in element]
    lines_sin_headers = lines[int(lines_header_CHROM[0])+1:]

    # Archivo para guardar las posiciones multialélicas
    multiallelic_vcf = open("{}.multiallelic.snp.vcf".format(prefix), "w+")
    multiallelic_positions = []

    # Gaurdaremos las filas multialélicas en el archivo .multiallelic.snp.vcf, y en el .remade mantendremos aquellas multialélicas que sean de SNPs (eg. G -> C,A) pero divididas en 2 filas para que las coja
    # el Minos. 
    with open("{}.remade.snp.vcf".format(prefix), "w+") as output_file:
        output_file.write("".join(headers))

        for line in lines_sin_headers: # Para cada posición
            tokens = line.split("\t")
            
            ALT = tokens[4]
            if "," not in ALT: # Si no es multialélica, la guardamos directamente
                output_file.write(line)
            
            else:
                multiallelic_positions.append(line.strip())
                # Primero comprobamos si alguna de las variantes posibles en la posicion no es un SNP
                ALT_list = ALT.split(",")
                lengths = [len(element) for element in ALT_list] # Saco las longitudes de cada variante posible
                lengths_greater_than_1 = [i for i in lengths if i > 1] # Me quedo con las que tienen len > 1 si las hubiera para comprobar si hay
                
                
                if len(lengths_greater_than_1) > 0: # si hay alguna variante que no sea un SNP (chunks), la descartamos
                    continue
                
                else: 

                    AS_SB_TABLE_info_list = tokens[7].split(";")[0][12:].split("|")
                    AS_SB_TABLE_info_list_list_to_write = []
                    for i in range(1,len(AS_SB_TABLE_info_list)):
                        AS_SB_TABLE_info_list_list_to_write.append([AS_SB_TABLE_info_list[0],AS_SB_TABLE_info_list[i]])
                    
                    DP_to_write = tokens[7].split(";")[1]
                    ECNT_to_write = tokens[7].split(";")[2]

                    MBQ_list = tokens[7].split(";")[3][4:].split(",")
                    MBQ_list_to_write = []
                    for i in range(1,len(MBQ_list)):
                        MBQ_list_to_write.append([MBQ_list[0],MBQ_list[i]])
                    
                    MFRL_list = tokens[7].split(";")[4][5:].split(",")
                    MFRL_list_to_write = []
                    for i in range(1,len(MFRL_list)):
                        MFRL_list_to_write.append([MFRL_list[0],MFRL_list[i]])
                    
                    MMQ_list = tokens[7].split(";")[5][4:].split(",")
                    MMQ_list_to_write = []
                    for i in range(1,len(MMQ_list)):
                        MMQ_list_to_write.append([MMQ_list[0],MMQ_list[i]])

                    MPOS_list_to_write = tokens[7].split(";")[6][5:].split(",")
                    POPAF_list_to_write = tokens[7].split(";")[7][6:].split(",")
                    TLOD_list_to_write = tokens[7].split(";")[8][5:].split(",")

                    GT_to_write = "0/1"

                    AD_list = tokens[9].split(":")[1].split(",")
                    AD_list_to_write = []
                    for i in range(1,len(AD_list)):
                        AD_list_to_write.append([AD_list[0],AD_list[i]])
                    
                    AF_list_to_write = tokens[9].split(":")[2].split(",")
                    DP_to_write = tokens[9].split(":")[3]

                    F1R2_list = tokens[9].split(":")[4].split(",")
                    F1R2_list_to_write = []
                    for i in range(1,len(F1R2_list)):
                        F1R2_list_to_write.append([F1R2_list[0],F1R2_list[i]])

                    F2R1_list = tokens[9].split(":")[5].split(",")
                    F2R1_list_to_write = []
                    for i in range(1,len(F2R1_list)):
                        F2R1_list_to_write.append([F2R1_list[0],F2R1_list[i]])

                    FAD_list = tokens[9].split(":")[6].split(",")
                    FAD_list_to_write = []
                    for i in range(1, len(FAD_list)):
                        FAD_list_to_write.append([FAD_list[0],FAD_list[i]])
                    
                    SB_to_write = tokens[9].split(":")[7]

                    for i in range(0,len(ALT_list)):
                        output_file.write('''{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tAS_SB_TABLE={7};DP={8};{9};MBQ={10};MFRL={11};MMQ={12};MPOS={13};POPAF={14};TLOD={15}\t{16}\t{17}:{18}:{19}:{20}:{21}:{22}:{23}:{24}'''.format(
                            tokens[0],tokens[1],tokens[2],tokens[3],ALT_list[i],tokens[5],tokens[6],("|").join(AS_SB_TABLE_info_list_list_to_write[i]), DP_to_write,
                            ECNT_to_write,(",").join(MBQ_list_to_write[i]),(",").join(MFRL_list_to_write[i]),(",").join(MMQ_list_to_write[i]),MPOS_list_to_write[i],
                            POPAF_list_to_write[i],TLOD_list_to_write[i], tokens[8], GT_to_write, (",").join(AD_list_to_write[i]),AF_list_to_write[i],DP_to_write,(",").join(F1R2_list_to_write[i]),
                            (",").join(F2R1_list_to_write[i]),(",").join(FAD_list_to_write[i]),SB_to_write
                        ))           
    
    # Guardamos las multialélicas
    multiallelic_vcf.write("".join(headers))
    multiallelic_vcf.write("\n".join(multiallelic_positions))
    multiallelic_vcf.close()
    UpdateHistory("Multiallelic filtering and trialelic readjudication for Mutect2 output", "custom", prefix)
    return 0

def callWT(prefix):
    '''Creates a Wild Type (WT) position file along with their read depth'''
    import vcf
    from .History import UpdateHistory
    vcf_reader = vcf.Reader(open("{}.gvcf".format(prefix), 'r'))

    wt = {}
    for record in vcf_reader:
        position, ref, alt  = record.POS, record.REF, record.ALT
        depth, genotype = record.samples[0]["DP"], record.samples[0]["GT"]
        alleles_depth = record.samples[0]["AD"]

	#check that there is no variant above 5% in the position
        if genotype == "0/0" and float(alleles_depth[0]) > (float(depth * 0.90)):
            wt[position] = depth

	#if genotype is not 0/0, again make sure there are no variants above 5%
	#frequency, and then report the correct depth for the WT allele
        elif len(alt) > 1:
            frequencies = record.samples[0]["AF"]

            if all([True if freq <= 0.10 else False for freq in frequencies]):
                wt[position] = alleles_depth[0]

    # save WT positions
    with open("{}.wt".format(prefix), 'w') as WT_file:
        WT_file.write("Position\tWT_Depth\n")
        for pos, d in wt.items():
            if d >= 6:
                WT_file.write(str(pos) + "\t" + str(d) + "\n")
    WT_file.close()
    UpdateHistory("WT obtained", "custom", prefix)


def Minos(reference,minos, snpEff, prefix, single_end):
    '''Perform final variant call adjudication with Minos'''
    import sys, os
    import subprocess as sp
    import shutil
    from subprocess import call
    from .History import UpdateHistory
    
    #first we get rid of multiallelic positions, as they're stored in another
    #file and Minos crashes if the input VCF has them
    
    if not single_end:
        cmd = ["singularity", "exec", "--bind", "/"+os.getcwd().split("/")[1],
            minos, "minos", "adjudicate", "--force",
            "--reads", "{}.P1.filtered.fastq.gz".format(prefix),
            "--reads", "{}.P2.filtered.fastq.gz".format(prefix),
            "{}_minos".format(prefix), reference.split("/")[-1],
            "{}.remade.snp.vcf".format(prefix), "{}.parsed.vcf".format(prefix)]
    else:
         cmd = ["singularity", "exec", "--bind", "/"+os.getcwd().split("/")[1],
            minos, "minos", "adjudicate", "--force",
            "--reads", "{}.filtered.fastq.gz".format(prefix),
            "{}_minos".format(prefix), reference.split("/")[-1],
            "{}.remade.snp.vcf".format(prefix), "{}.parsed.vcf".format(prefix)]       
           
    stat = call(cmd)
    if stat != 0:
        sys.exit("Pipeline stopped at Minos variant adjudication for sample {}!\n".format(prefix))
    UpdateHistory(cmd, "minos", prefix)

    # We extract the Minos output from each folder
    sp.run("cd {0}_minos; grep -v '0/0' final.vcf > {1}.final_sin_wt.vcf; mv {2}.final_sin_wt.vcf ..; cd ..".format(prefix,prefix,prefix), shell=True)
    #sp.run("find . -print | grep '.final_sin_wt.vcf' | xargs -I {} -P 20 cp {} .", shell=True)
    UpdateHistory("WT filtered from Minos output", "custom", prefix)
    # We annotate the Minos output with SnpEff

    cmd = ["java", "-jar", snpEff, "-v", "-hgvs1LetterAa", "-noStats", "-no-downstream", "-no-intergenic", "-no-intron", "-no-upstream", "-noLof",
           "MTB_ancestor", "{}.final_sin_wt.vcf".format(prefix)]

    with open("{}.final_sin_wt.annotSnpEff.vcf".format(prefix), "w+") as outfh:
        stat = call(cmd, stdout=outfh)
    UpdateHistory("Anotado con snpEff: java -jar snpEff -v -hgvs1LetterAa -noStats -no-downstream -no-intergenic -no-intron -no-upstream -noLof MTB_ancestor {}.final_sin_wt.vcf".format(prefix),"custom",prefix)


def minos_raw_vcf_to_tab(prefix):
    # Convierte el output del minos a tab quitando posiciones Cov_total < 3 y VarFreq < 0.05. Para las variantes bialélicas con una freq de ALT < 90, usa el diccionario IUPAC para la columna Cons y añade la variante original a la última columna
    # Para bialélicas con una freq de ALT > 90, mantiene la información de la posición y añade la variante a la última columna
    # Para trialélicas, si alguna está a maś del 90%, la usa como variante en dicha línea
    # Para trialélicas, si no hay ninguna fijada Y hay una con más DP que las demás, calcula la freq para todas las variantes propuestas
    # y en Cons añade una "?", pero mantiene las variantes en la última columna ordenadas de igual manera que las frecuencias para que podamos
    # saber la freq de cada variante.
    # Para trialélicas que tienen 2 variantes con DP máxima (ej 2,6,6), busca su frecuencia en los raw de varscan y mutect2 y las añade en el mismo
    # orden que las variantes propuestas añadidas en la última columna. En la columna Cons también pondrá un "?"
    #from .History import UpdateHistory
    import subprocess as sp
    iupac = {"M": ["A", "C"],
        "R": ["A", "G"],
        "W": ["A", "T"],
        "S": ["C", "G"],
        "Y": ["C", "T"],
        "K": ["G", "T"]}
    
    def find_key_with_elements(d, elements):
        '''Función para sacar la clave que contenga los elementos que le indiquemos'''
        required_elements = set(elements)
        for key, value in d.items():
            if required_elements.issubset(set(value)):
                return key
        return None

    def get_freqs_from_varscan_and_mutect(pos_minos):
        '''Esta función obtendrá la frecuencia para cada posición consensuada por el Minos a partir de los output de VarScan y Mutect2. Si la posición
           está en ambos callers, se hará la media. Si no, se usará la que esté disponible'''
        
        variantes = cons.split(",")
        lista_keys = [ref+i for i in variantes]
        dic_variants = {key: [] for key in lista_keys}

        # Ahora vamos a buscar las freqs en el archivo del varscan y del mutect
        with open("{}.parsed.vcf".format(prefix), "r+") as input_varscan:
            lines_varscan = input_varscan.readlines()
    
            for line_varscan in lines_varscan:

                try:
                    if "MTB_anc\t{}\t".format(pos_minos) in line_varscan: # Si está la posición, añadimos la freq a la clave del diccionario correspondiente
                        ref_varscan, cons_varscan, freq_varscan = itemgetter(3,4,9)(line_varscan.split("\t"))
                        freq_variant_varscan = freq_varscan.split(":")[-1] # nos quedamos con la frecuencia de la variante
                        dic_variants[ref_varscan+cons_varscan].append(float(freq_variant_varscan)) # Guardamos la frecuencia de la variante en su clave
                except KeyError: # Puede ocurrir que alguna posición tenga 2 posibles variantes en VarScan. Sin embargo, si alguna no aparece en la columna ALT
                                 # de Minos, no se va a tener en cuenta
                    continue

                    # Ahora hacemos lo mismo con Mutect2
        with open("{}.remade.snp.vcf".format(prefix), "r+") as input_mutect2:
            lines_mutect2 = input_mutect2.readlines()

            for line_mutect2 in lines_mutect2:
                try:
                    if "MTB_anc\t{}\t.".format(position) in line_mutect2:
                        ref_mutect, cons_mutect, info_mutect = itemgetter(3,4,9)(line_mutect2.split("\t"))
                        freq_variant_mutect = info_mutect.split(":")[2] # Guardamos la frecuencia
                        dic_variants[ref_mutect+cons_mutect].append(float(freq_variant_mutect))
                except KeyError:
                    continue
                    
        from statistics import mean

        # Guardamos la media de las frecuencias de cada variante extraídas de VarScan y Mutect2 y se añaden a la columna de VarFreq ordenadas de mayor a menor
        for keys,values in dic_variants.items():
            if len(values) == 1:
                dic_variants[keys].append(0)

        freqs_to_write = [round(mean(values)*100,2) for keys,values in dic_variants.items()]
        if len(freqs_to_write) == 1:
            return freqs_to_write[0]
                    
        elif len(freqs_to_write) > 1:
            return freqs_to_write


    def get_depths_from_varscan_and_mutect(pos_minos):
        '''Esta función obtendrá la profundidad para cada posición consensuada por el Minos a partir de los output de VarScan y Mutect2. Si la posición
           está en ambos callers, se hará la media. Si no, se usará la que esté disponible'''
        
        variantes = cons.split(",")
        lista_keys = [ref+i for i in variantes]
        dic_variants = {key: [] for key in lista_keys}

        # Ahora vamos a buscar las freqs en el archivo del varscan y del mutect
        with open("{}.parsed.vcf".format(prefix), "r+") as input_varscan:
            lines_varscan = input_varscan.readlines()
    
            for line_varscan in lines_varscan:

                try:
                    if "MTB_anc\t{}\t".format(pos_minos) in line_varscan: # Si está la posición, añadimos la freq a la clave del diccionario correspondiente
                        ref_varscan, cons_varscan, freq_varscan = itemgetter(3,4,9)(line_varscan.split("\t"))
                        depth_variant_varscan = freq_varscan.split(":")[1] # nos quedamos con la depth de la variante
                        dic_variants[ref_varscan+cons_varscan].append(float(depth_variant_varscan)) # Guardamos la frecuencia de la variante en su clave
                except KeyError: # Puede ocurrir que alguna posición tenga 2 posibles variantes en VarScan. Sin embargo, si alguna no aparece en la columna ALT
                                 # de Minos, no se va a tener en cuenta
                    continue

                    # Ahora hacemos lo mismo con Mutect2
        with open("{}.remade.snp.vcf".format(prefix), "r+") as input_mutect2:
            lines_mutect2 = input_mutect2.readlines()

            for line_mutect2 in lines_mutect2:
                try:
                    if "MTB_anc\t{}\t.".format(position) in line_mutect2:
                        ref_mutect, cons_mutect, info_mutect = itemgetter(3,4,9)(line_mutect2.split("\t"))
                        depth_variant_mutect = info_mutect.split(":")[1].split(",")[1] # Guardamos la depth de la variante
                        dic_variants[ref_mutect+cons_mutect].append(float(depth_variant_mutect))
                except KeyError:
                    continue
                    
        from statistics import mean

        # Guardamos la media de las frecuencias de cada variante extraídas de VarScan y Mutect2 y se añaden a la columna de VarFreq ordenadas de mayor a menor
        for keys,values in dic_variants.items():
            if len(values) == 1:
                dic_variants[keys].append(0)

        depths_to_write = [int(mean(values)) for keys,values in dic_variants.items()]

        if len(depths_to_write) == 1:
            return depths_to_write[0]
                    
        elif len(depths_to_write) > 1:
            return depths_to_write



    def check_bigger_than_90(list1, val):

        for x in list1:
            if x >= val:
                return True
    
    def check_lower_than_90(list1, val):
     
    # traverse in the list
        for x in list1:
            if x < val:
                return True


    
    from operator import itemgetter

    with open("{}.final_sin_wt.vcf".format(prefix), "r+") as filtered_raw:
        lines = filtered_raw.readlines()

    with open("{}.minos.raw.tab".format(prefix), "w+") as raw_tab:
        raw_tab.write("#Chrom\tPosition\tRef\tCons\tVarFreq\tCov_allele\tVarAllele\n") # Escribimos los header del archivo tabla

        # Quitamos los headers del filtered_raw
        lines_noheader = [i for i in lines if "#" not in i]

        for line in lines_noheader: # Guardamos la info que nos interesa
            tokens = line.split("\t")
            chrom, position, ref, cons = itemgetter(0,1,3,4)(tokens)

            if "./.:.:0,0:.:0:0,0:0.0:0.0" in line:
                chrom = "MTB_anc*"

            varfreq = get_freqs_from_varscan_and_mutect(position)
            cov_total = get_depths_from_varscan_and_mutect(position)

            if len(cons.split(",")) == 1 and varfreq < 90: # Si solo hay 1 variante y no esta fijada, usamos la nomenclatura IUPAC y añadimos la variante a VarAllele
                code_iupac = find_key_with_elements(iupac, [ref, cons])
                raw_tab.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,ref,code_iupac,varfreq,cov_total,cons))
    
            elif len(cons.split(",")) == 1 and varfreq >= 90: # Si está fijada, la mantenemos
                raw_tab.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,ref,cons,varfreq,cov_total,cons))

            elif len(cons.split(",")) > 1 and check_bigger_than_90(varfreq,90): # Si hay más de una y alguna está al menos al 90%, mantenemos esa en Cons pero añadimos todas las opciones y freqs
                variantes = cons.split(",")
                # Sacamos el indice del valor de freq más alto entre las variantes posibles, que corresponderá a la variante con mayor freq
                max_index = varfreq.index(max(varfreq))
                min_index = varfreq.index(min(varfreq))
                variante_90perc = variantes[max_index] # Sacamos la variante > 90%
                variante_minoritaria = variantes[min_index] 
                string_varfreq = "{},{}".format(varfreq[max_index], varfreq[min_index])
                string_VarAllele = "{},{}".format(variante_90perc, variante_minoritaria)
                string_depth = "{},{}".format(cov_total[max_index], cov_total[min_index])
                
                # Ahora, añadimos la línea pero solo manteniendo la variante > 90 % en Cons
                raw_tab.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,ref,variante_90perc,string_varfreq,string_depth,string_VarAllele))
            
            # Si no hay ninguna variante fijada
            elif len(cons.split(",")) > 1 and check_lower_than_90(varfreq,90): 

                variantes = cons.split(",")

                # Sacamos el indice del valor de freq más alto entre las variantes posibles, que corresponderá a la variante con mayor freq
                max_index = varfreq.index(max(varfreq))
                min_index = varfreq.index(min(varfreq))
                variante_mayoritaria = variantes[max_index] # Sacamos la variante a mayor freq
                variante_minoritaria = variantes[min_index] 
                string_varfreq = "{},{}".format(varfreq[max_index], varfreq[min_index])
                string_VarAllele = "{},{}".format(variante_mayoritaria, variante_minoritaria)
                string_depth = "{},{}".format(cov_total[max_index], cov_total[min_index])

                raw_tab.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,ref,"?",string_varfreq,string_depth,string_VarAllele))
                                   
    filter_raw_minos("{}.minos.raw.tab".format(prefix),prefix)
    return 0

def filter_raw_minos(file_name,prefix):
    '''Función para filtrar el raw del minos en formato tabla manteniendo las posiciones con una depth > 3 y freq > 5'''
    from operator import itemgetter
    #from .History import UpdateHistory
    # Abrimos el archivo a filtrar
    with open(file_name, "r+") as input_file:
        lines_input = input_file.readlines()

    # Generamos el archivo ya filtrado

    posiciones_a_mantener = []

    with open("{}.filtered.minos.raw.tab".format(file_name.split(".")[0]), "w+") as output_file:
        # Mantenemos los headers
        headers = [line for line in lines_input if "#" in line]
        lines_noheader = [line for line in lines_input if "#" not in line]
        output_file.write(("").join(headers))

        # Ahora filtramos aquellas filas con un COV_TOTAL < 3

        for line in lines_noheader: # Iteramos desde la linea siguiente del header hasta el final del archivo
            depth, freq = itemgetter(5,4)(line.split("\t"))
            depth = depth.split(",")[0] # Si hubiera más de una profundidad, nos quedamos con la mayor, que es la que está primero
            freq = freq.split(",")[0] # Si hubiera más de una freq, nos quedamos con la mayor, que es la que está primero
            if int(depth) >= 3 and float(freq) >= 5: # Si la posicion tiene más de 3 lecturas y una frecuencia de al menos 0.05, la mantenemos
                posiciones_a_mantener.append(line)
            else:
                continue

        output_file.write(("").join(posiciones_a_mantener)) # Guardamos las posiciones a mantener
    #UpdateHistory("Minos output filtered keeping variants with depth > 3 and freq > 5", "custom", prefix)
    return 0

  
# EPI and density filtering
def filter_EPI(prefix):
    from .History import UpdateHistory
    '''Dejamos aquellas variantes del output del Minos que tengan al menos 20 de coverage y 90% de frecuencia'''
    
    with open("{}.filtered.minos.raw.tab".format(prefix), "r+") as input_file:
        lines = input_file.readlines()
    
    with open("{}.EPI.snp.nodensityfilter".format(prefix), "w+") as output:
        output.write(lines[0]) # Escribimos el header
        lines = lines[1:] # Quitamos el header

        for line in lines:
            tokens = line.split("\t")
            cons, coverage, var_freq, var_allele = tokens[3], tokens[5], tokens[4], tokens[6]

            if "?" in cons: # Si la variante es trialelica y no hay ninguna fijada, pasamos a la siguiente línea
                continue

            elif "," in var_allele and cons in ["A","T","C","G"]: # Si es trialélica, miramos si alguna está fijada y guardamos esa
                highest_coverage = coverage.split(",")[0] # Si hubiera más de una profundidad, nos quedamos con la mayor, que es la que está primero
                highest_var_freq = var_freq.split(",")[0] # Si hubiera más de una freq, nos quedamos con la mayor, que es la que está primero 
                main_variant = tokens[6].split(",")[0]

                if int(highest_coverage) >= 20 and float(highest_var_freq) >= 90:
                    output.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(tokens[0],tokens[1],tokens[2],main_variant,highest_var_freq,
                                                                            highest_coverage,main_variant))
                else:
                    continue

            elif int(coverage) >= 20 and float(var_freq) >= 90:
                output.write(line)
    
    return 0


def densityfilter(prefix, window, dens):
    '''Filter those SNPs that are in high density areas of SNPs'''
    from .History import UpdateHistory
    import pandas as pd
    pos_to_remove = []

    EPI_tab_file = pd.read_csv("{}.EPI.snp.nodensityfilter".format(prefix), sep = '\t')

    # header to inherit
    with open("{}.EPI.snp.nodensityfilter".format(prefix), 'r') as f:
        header = [li for li in f if li.startswith('#')]

    # for each POS
    # how many other POS are less than 'window' pb away from them
    # Is this number greater than dens?
    for p in EPI_tab_file['Position']:
        aux = list(filter(lambda x:
                          x >= (p-window) and x <= (p+window),
                          EPI_tab_file['Position']))
        if len(aux) >= dens:
            pos_to_remove.extend(aux)

    # remove duplicates and sort
    pos_to_remove = list(set(pos_to_remove))
    pos_to_remove.sort()

    # save positions removed by density
    with open("{}.dens_removed_snps.tab".format(prefix), 'w') as rm_file:
        rm_file.write('\t'.join(header))
        EPI_tab_file.loc[EPI_tab_file['Position'].isin(pos_to_remove)].to_csv("{}.dens_removed_snps.tab".format(prefix),mode='a', sep="\t", index=False)


    # save filtered VCF
    final_file = open("{}.EPI.snp.final".format(prefix), 'w')
        #final_file.write('\t'.join(header))

    EPI_tab_file.loc[~EPI_tab_file['Position'].isin(pos_to_remove)
                 ].to_csv("{}.EPI.snp.final".format(prefix),
                          mode='a', sep="\t", index=False)
    final_file.close()

    print(len(pos_to_remove), " SNPs removed by density filter: "
          "window = ", window, ", density = ", dens)

    #UpdateHistory("Density filter applied to {}.EPI.snp.vcf file"
    #              " with parameters window = {} and density = "
    #              "{}".format(prefix, window, dens),
    #              "custom", prefix)


def LoadAnnotation(anno_file):
    '''Load the annotation file to store each position to
    be discarded in a list

    This function is meant to load by default
    the H37Rv annotation under
    /data/Databases/MTB_annotation/H37Rv.annotation_new.tsv
    '''
    from .History import UpdateHistory
    discard_positions = []
    with open(anno_file) as infile:
        next(infile)
        for line in infile:
            line = line.strip().split("\t")

            if line[-1] == "DISCARD":
                start, end = int(line[1]), int(line[2])
                for pos in range(start, end + 1):
                    discard_positions.append(str(pos))

    return discard_positions


def FilterSnps(prefix):
    '''Filter SNP files according to the TSV file
    passed. This TSV file is an annotation file of the
    H37Rv genome, with the last column indicating if SNPs in 
    this region must be kept or not'''
    
    import os
    from .History import UpdateHistory
    annotation = LoadAnnotation("/data/ThePipeline_v2/data/H37Rv.annotation_new.tsv".format(
                    os.path.split(
                        os.path.dirname(
                            os.path.abspath(__file__)))[0]))

    with open("{}.EPI.snp.final.annoF".format(prefix), "w") as outfile:
        with open("{}.EPI.snp.final".format(prefix)) as infile:
            header = infile.readline()
            outfile.write(header)
            for line in infile:
                pos = line.strip().split("\t")[1]
                if pos not in annotation:
                        outfile.write(line)

    return 0

def correct_line(content):

    iupac = {"M": ["A", "C"],
        "R": ["A", "G"],
        "W": ["A", "T"],
        "S": ["C", "G"],
        "Y": ["C", "T"],
        "K": ["G", "T"]}

    amb = iupac[content[3]]

    if (content[2] in amb) and (content[-1] in amb):
        return(content)
    else:
        amb.remove(content[2])
        content[-1] = amb[0]
        return(content)
	
	
def get_DR(prefix):
    '''Perform annotation of DR.snp files and lowcov positions'''
    import os
    from .AnnotatePos import annotateSingle, annotateDouble, annotateTriple

    annotation = "/data/ThePipeline_v2/data/annotation_H37Rv.csv"
    
    ancestor = open("/data/Databases/MTB_ancestor/MTB_ancestor_reference.fasta").read().replace("\n","").replace(">MTB_anc","")
    
    code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"!", "TAG":"!",
    "TGT":"C", "TGC":"C", "TGA":"!", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    genes, names, output, codons = {}, {}, [], []
    
    #first we read the annotation
    with open(annotation) as inAnnot:
        for line in inAnnot:
            line = line.strip().split(",")
            genes[line[0]] = [int(line[-3]),int(line[-2]),line[-1]]
            names[line[0]] = line[1]

    import sys
    #print(genes)

    print("\033[92m\nNow annotating variants from sample {}...\n\033[00m".format(prefix))
    #now we load the DR-snp file and annotate each individual position
    
    with open("{}.filtered.minos.raw.tab".format(prefix)) as infile:
        next(infile)

        for line in infile:
            line = line.strip().split("\t")

            if "," in line[4]: # Modificar esto para que despliegue las líneas con trialélicas en 2 líneas separadas

                freqs = line[4].split(",")
                coverage = line[5].split(",")
                variants = line[6].split(",")
                line_split = [[line[0], line[1], line[2], variants[0], freqs[0], coverage[0], variants[0]],[line[0], line[1], line[2], variants[1], freqs[1], coverage[1], variants[1]]]

                for line_ in line_split:
                # Solo mantenemos las posiciones con una freq mayor de 5% y una profundidad mayor de 6
                    if float(line_[4]) < 5 and float(line_[5]) < 6:
                        continue
                    if line_[3] in ["M", "R", "W", "S", "Y", "K"]:
                        line_ = correct_line(line_)      

                    if len(line_[-1]) == 1:
                        snpid = line_[1]+line_[2]+line_[-1] # POS + REF + ALT
                        change = annotateSingle(snpid, genes, code, ancestor)
                    # print(line)
                        line_.extend(change.strip().split(","))
                        line_.append("triallelic")
                    # print(line)
                        output.append(line_)
            
            else:
                if float(line[4]) < 5 and float(line[5]) < 6:
                    continue
                if line[3] in ["M", "R", "W", "S", "Y", "K"]:
                    line = correct_line(line)      

                if len(line[-1]) == 1:
                    snpid = line[1]+line[2]+line[-1] # POS + REF + ALT
                    change = annotateSingle(snpid, genes, code, ancestor)
                    # print(line)
                    line.extend(change.strip().split(","))
                    # print(line)
                    output.append(line)



    for item in output:
        if "intergenic" not in item:
            codons.append(item[7]+"_"+item[8][1:-1]) # Se guarda la posicion del codon
        else:
            codons.append(item[7])
    
    
    #here we check if there is more than one mutation in the same codon
    double_muts, triple_muts = [], []
    for i in range(0, len(codons)):
        try:
            if codons[i] == codons[i+1]: #we have a multiple mutation

                if codons[i] == codons[i+2]: #this is a triple mutation
                    if "IG" not in codons[i]:
                        triple_muts.append(codons[i])
                        
                else: #it is a double then
                    if "IG" not in codons[i]:
                        double_muts.append(codons[i])

        except IndexError:
            pass
    
    double_muts, triple_muts = list(set(double_muts)), list(set(triple_muts))
    for item in triple_muts:
        if item in double_muts:
            double_muts.remove(item)

    final_output, computed = [], {}
    for i in range(0, len(output)): #now let's correct the preliminary annotation accordingly
        this_codon = output[i][7]+"_"+output[i][8][1:-1]

        if this_codon in double_muts:
            if this_codon not in computed.keys():

                first_mut = output[i][1]+output[i][2]+output[i][6]
                second_mut = output[i+1][1]+output[i+1][2]+output[i+1][6]
                double_change = annotateDouble(first_mut + "+" + second_mut, genes, code, ancestor)

                computed[this_codon] = double_change
            else:
                double_change = computed[this_codon]
                
            new_line = output[i]
            new_line[7:] = double_change.strip().split(",")
            if output[i][1] == output[i+1][1] or output[i][1] == output[i-1][1]:
                new_line.append("TRIALLELIC POSITION")
            else:
                new_line.append("Double mutation")
            final_output.append(new_line)
            
        elif this_codon in triple_muts:
            if this_codon not in computed.keys():
                first_mut = output[i][1]+output[i][2]+output[i][6]
                second_mut = output[i+1][1]+output[i+1][2]+output[i+1][6]
                third_mut = output[i+2][1]+output[i+2][2]+output[i+2][6]

                triple_change = annotateTriple(first_mut + "+" + second_mut + "+" + third_mut, genes, code, ancestor)
                computed[this_codon] = triple_change
            else:
                triple_change = computed[this_codon]
                
            new_line = output[i]
            new_line[7:] = triple_change.strip().split(",")
            if output[i+1][1] == output[i+2][1] or output[i][1] == output[i+1][1] or output[i][1] == output[i-1][1]:
                new_line.append("Double mutation (TRIALLELIC)")
            else:
                new_line.append("Triple mutation")
            final_output.append(new_line)
        else:
            final_output.append(output[i])

    #finally we save the resulting annotated DR.snp.final file
    with open("{}.DR.snp.final".format(prefix), "w") as outfile:
        outfile.write("\t".join(['Chrom', 'Position', 'Ref', 'Cons', 'VarFreq', 'Cov_total', 'VarAllele','Gene','Change','Ann','CodonWT','CodonVAR','Comments\n']))
        for line in final_output:
            outfile.write("\t".join(line)+"\n")
    
    low_cov = {}
    with open("{}.lowcov.tsv".format(prefix)) as infile:
        next(infile)
        for line in infile:
            line = line.strip().split("\t")
            this_gene = ""
            for gene, coords in genes.items():
                if coords[0] <= int(line[0]) <= coords[1]:
                    this_gene = gene
            low_cov[line[0]] = this_gene
            
    with open("{}.lowcov".format(prefix), "w") as outfile:
        outfile.write("Position\tLocus\tGene\n")
        for position, gene in low_cov.items():
            outfile.write("\t".join([position, gene, names[gene]+"\n"]))

def CorrectDR(prefix):
    '''Función para corregir la anotación de las posiciones trialélicas'''

    from .AnnotatePos import annotateSingle, annotateDouble, annotateTriple

    annotation = "/data/ThePipeline_v2/data/annotation_H37Rv.csv"
    
    ancestor = open("/data/Databases/MTB_ancestor/MTB_ancestor_reference.fasta").read().replace("\n","").replace(">MTB_anc","")
    
    code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"!", "TAG":"!",
    "TGT":"C", "TGC":"C", "TGA":"!", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    genes, names, output, codons = {}, {}, [], []
    
    #first we read the annotation
    with open(annotation) as inAnnot:
        for line in inAnnot:
            line = line.strip().split(",")
            genes[line[0]] = [int(line[-3]),int(line[-2]),line[-1]]
            names[line[0]] = line[1]

    with open("{}.DR.snp.final".format(prefix), "r+") as DR_file:
        lines = DR_file.readlines()

    with open("{}.DR.snp.final".format(prefix), "w+") as output_DR_file:
        output_DR_file.write(lines[0])

        to_skip = []

        for i in range(1, len(lines)):
            line = lines[i].strip().split("\t")

            if line[1] in to_skip:
                continue
            if "TRIALLELIC POSITION" in line: # Si la trialélica no está en una doble mutación con otra posición del mismo codón, la reanotamos de manera independiente para cada alelo
                new_line = line[0:7]
                snpid = new_line[1]+new_line[2]+new_line[-1] # POS + REF + ALT
                change = annotateSingle(snpid, genes, code, ancestor)
                new_line.extend(change.strip().split(","))
                new_line.append("TRIALLELIC POSITION")
                output_DR_file.write("\t".join(new_line) + "\n")
            
            elif "Double mutation (TRIALLELIC)" in line:

                first_mut = line[1]+line[2]+line[6]
                second_mut = lines[i+1].split("\t")[1]+lines[i+1].split("\t")[2]+lines[i+1].split("\t")[6]
                first_double_change = annotateDouble(first_mut + "+" + second_mut, genes, code, ancestor)
                second_mut = lines[i+2].split("\t")[1]+lines[i+2].split("\t")[2]+lines[i+2].split("\t")[6]
                second_double_change = annotateDouble(first_mut + "+" + second_mut, genes, code, ancestor)

                new_line_first = line[0:7]
                new_line_first.extend(first_double_change.strip().split(","))
                new_line_first.append("Double mutation (TRIALLELIC)")
                new_line_second = line[0:7]
                new_line_second.extend(second_double_change.strip().split(","))
                new_line_second.append("Double mutation (TRIALLELIC)")

                new_line_third = lines[i+1].split("\t")[0:7]
                new_line_third.extend(first_double_change.strip().split(","))
                new_line_third.append("Double mutation (TRIALLELIC)")

                new_line_fourth = lines[i+2].split("\t")[0:7]
                new_line_fourth.extend(second_double_change.strip().split(","))
                new_line_fourth.append("Double mutation (TRIALLELIC)")

                to_skip.append(line[1])
                to_skip.append(lines[i+1].split("\t")[1])

                output_DR_file.write("\t".join(new_line_first) + "\n" + "\t".join(new_line_third) + "\n" + "\t".join(new_line_second) + "\n" + "\t".join(new_line_fourth) + "\n")
            
            else:
                output_DR_file.write("\t".join(line)+"\n")


def annotateVCF(snpEff, prefix):
    '''Annotate the VCF files using SnpEff'''
    import subprocess as sp
    import sys
    from os.path import exists
    from .History import UpdateHistory

    cmd_snpEff = ["java", "-jar", snpEff, "-v", "-hgvs1LetterAa",
                  "-noStats", "-no-downstream", "-no-intergenic",
                  "-no-intron", "-no-upstream", "-noLof",
                  "MTB_ancestor", "{}.EPI.snp.vcf".format(prefix)]

    with open("{}.annot".format(prefix), "w") as outfh:
        stat = sp.call(cmd_snpEff, stdout=outfh)

    # rename so the annotated VCF will be the original
    sp.run("mv {}.annot {}.EPI.snp.vcf".format(prefix, prefix),
               shell=True, capture_output=True)

    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}:"
                 " Pipeline stopped while annotating "
                 "with SnpEff!\n\033[0m".format(prefix))

    UpdateHistory(cmd_snpEff, "snpeff", prefix)
    
    
def Calling(args):
    from .Repository import Programs, Data
    import sys
    import os
    import subprocess as sp
    import shutil
    # Check package dependencies
    try:
        import vcf
    except ModuleNotFoundError:
        sys.exit("\033[91mPIPELINE ERROR: The VCF module (PyVCF) is not installed "
                 "and is necessary for the 'Calling' process\033[0m")
    try:
        import pandas
    except ModuleNotFoundError:
        sys.exit("\033[91mPIPELINE ERROR: The PANDAS module is not installed "
                 "and is necessary for the 'Calling' process\033[0m")

    programs = Programs()
    data = Data()

    varscan = programs["varscan"]
    minos = programs["minos"]
    gatk = programs["gatk"]
    snpEff = programs["snpEff"]
    samtools = programs["samtools"]
    reference = args.reference
    genomeCoverageBed = programs["genomeCoverageBed"]

    # Check that mapping files end with .bam/.cram
    check_extension(args.ext)

    if not reference:
        reference = data["reference"]
        
    # SNP calling with VarScan2
    print("\033[92m\nRunning VarScan on {}\n\033[00m".format(args.prefix))
    VarScan(reference, args.prefix, varscan, samtools, args.ext)

    # SNP calling with Mutect2
    print("\033[92m\nRunning Mutect2 on {}\n\033[00m".format(args.prefix))
    Mutect2(reference, args.prefix, gatk, samtools, genomeCoverageBed,
            args.ext, args.threads, args.min_depth, args.min_qual,
            args.filter_freq)

    # Multiallelic filter from Mutect2 snp output
    filter_multiallelic_from_mutect2_snp(args.prefix)

    #let's call the wild type positions from the gVCF file
    print("\033[92m\nCalling WT for {}\n\033[00m".format(args.prefix))
    callWT(args.prefix)

    # Minos calling

    #we now copy the reference sequence for Minos analysis to avoid
    #unnecessary singularity path bindings, and run the program inside its container
    os.popen('cp ' + reference + ' .')
    print("\033[92m\nPerforming variant call adjudication with Minos for sample {}...\n\033[00m".format(args.prefix))
    Minos(reference,minos,snpEff,args.prefix, args.single_end)
    minos_raw_vcf_to_tab(args.prefix)
    # EPI filtering
    print("\033[92m\nObtaining {}.EPI.snp.final\n\033[00m".format(args.prefix))
    filter_EPI(args.prefix)
        
    # Apply density filter
    if not args.skip_dens_filt:  # by density
        print("\033[92m\nFiltering {}.EPI.snp.final by density\n\033[00m".format(args.prefix))
        densityfilter(args.prefix, args.window, args.density)

    # annoF filtering
    print("\033[92m\nObtaining {}.EPI.snp.final.annoF\n\033[00m".format(args.prefix))
    FilterSnps(args.prefix)

    # DR extraction and correction
    print("\033[92m\nObtaining {}.DR.snp.final\n\033[00m".format(args.prefix, args.prefix))
    get_DR(args.prefix)
    CorrectDR(args.prefix)
    
    print("\033[92m\nRemoving and renaming files\n\033[00m")
    #a final cleanup and renaming of intermediate and/or useless files
    os.remove("{}.EPI.snp.nodensityfilter".format(args.prefix))
    os.rename("{}.filtered.minos.raw.tab".format(args.prefix), "{}.snp.minos".format(args.prefix))
    os.remove("{}.final_sin_wt.vcf".format(args.prefix))
    os.remove("{}.lowcov.tsv".format(args.prefix))
    os.remove("{}.minos.raw.tab".format(args.prefix))
    os.remove("{}.parsed.vcf".format(args.prefix))
    os.rename("{}.remade.snp.vcf".format(args.prefix),"{}.snp.mutect".format(args.prefix))
    #os.remove("{}.res".format(args.prefix))
    #os.rename("{}.res.corrected".format(args.prefix),"{}.res".format(args.prefix))
    os.rename("{}.snp".format(args.prefix),"{}.snp.varscan".format(args.prefix))
    os.remove("{}.snp.vcf".format(args.prefix))
    os.remove("{}.vcf".format(args.prefix))
    os.remove("{}.vcf.stats".format(args.prefix))

    if not args.keep_gvcf:
        os.remove("{}.gvcf".format(args.prefix))
        os.remove("{}.gvcf.stats".format(args.prefix))
    else:
        print("\033[92m\nCompressing {}.gvcf\n\033[00m".format(args.prefix))
        sp.run("gzip {}.gvcf".format(args.prefix), shell=True)

    if not args.keep_bam:
        # If .cram file not present, convert bam into cram and then remove bam
        if "{}.cram".format(args.prefix) not in os.listdir():
            cmd_toCRAM = [samtools, "view", "--threads", args.threads,
            "-T", reference,
            "-C", "-o", "{}.cram".format(args.prefix),
            "{}{}".format(args.prefix, args.ext)]
            sp.call(cmd_toCRAM)
            os.remove("{}.sort.bam".format(args.prefix))
            os.remove("{}.sort.bam.bai".format(args.prefix))
        else:
            os.remove("{}.sort.bam".format(args.prefix))
            os.remove("{}.sort.bam.bai".format(args.prefix))

    shutil.rmtree("{}_minos".format(args.prefix))
    
    return 0

