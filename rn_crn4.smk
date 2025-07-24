'''
This pipeline is focused in ring nucleases and transcription factors in type III CRISPR-Cas loci. 
Built on top of new_effectors.smk from a previous project.

Remember to set file paths below to fit your local machine after downloading databases as instructed in the README.md file.
'''

project = "run1" # name of the project. This will be used to create a folder for the results.

thread_hogger = 50 #number of threads dedicated to a single thread-hogging rule. Meant for rules that do not rely on wildcards.
thread_small = 5 #for wildcard tasks that benefit from multithreading
thread_ultrasmall = 1 #for single-thread wildcard based rules

ncbi_email = "your_email_here" # your email here, used for NCBI queries

local_protein_blast_db = "" #path to local protein blast database. For example, "/mnt/shared/apps/databases/ncbi/nr". Check README.

base_path = "your_base_path" + "/" + project # replace with your base path where the results will be stored
program_root = "your_program_path" # replace with where this pipeline is located on your local machine

cas10_cluster_threshold = 0.9 #initial Cas10 clustering threshold
crispr_locus_interference_cutoff = 0 #cutoff for CRISPR loci interference completeness. Loci with less than this percentage of interference genes present are discarded

protein_clustering = str(config["protein_clustering"]) # "True" or "False". If True, proteins are clustered before further analysis.
getGenomesBy = str(config["getGenomesBy"]) # "local" or "remote". If "local", genomes are fetched from a local folder. If "remote", genomes are fetched from NCBI using the NCBI dataset API.
cas10_anchor = config["cas10_anchor"] # "True" or "False". If True, Cas10 is used as an anchor for the pipeline, meaning that the pipeline will only process loci with Cas10 present. If False, all loci are processed regardless of Cas10 presence.
catyper_hmm_evalue = "1e-10" # e-value threshold for cATyper hmmsearch

hmm_msa_folder = "data/known_effectors/profiles" #folder names within this folder are used for creating effector dictionary
hmm_database_folder = "data/known_effectors/concatenated_profiles" #contains the concatenated and hmmpressed hmm profiles
hmm_database_file = "all.hmm" #filename for the concatenated hmmpressed hmm profiles

cas10_db = "data/cas10/all_cas10s.hmm" # hmm profile for Cas10 proteins
modified_cas10_hd = "data/cas10/HD_HMM.msa" #modified cas10 hmm profile for hmmsearch

TM_path = "data/TMHMM"
validated_effectors_hmm_db = "data/validated_new_effectors/concatenated_profiles/all.hmm" # hmm profile for validated new effectors
validated_effectors_folder = "/data/validated_new_effectors" # folder containing validated new effectors hmm profiles

ring_nuclease_db = "data/rns/all.hmm" # hmm profile for ring nucleases
ring_nuclease_folder = "data/rns" # folder containing ring nuclease hmm profiles

#public databases (local)
cogs_db = "" # path to your local COG database. E.g. /home/vhoikkal/scratch/private/databases/cogs/COG_KOG
pfam_db = "" # path to your local Pfam database. E.g. /home/vhoikkal/scratch/private/databases/pfam/Pfam-A.hmm
pdb30_db = "" # path to your local PDB30 database. E.g. /home/vhoikkal/scratch/private/databases/pdb30

millard_fa = "" # path to the Millard phage genomes file, e.g. "/mnt/shared/scratch/vhoikkal/private/databases/millard_phages/4May2024_genomes.fa"
millard_proteins = "" # path to the Millard phage proteins file, e.g. "/mnt/shared/scratch/vhoikkal/private/databases/millard_phages/4May2024_vConTACT2_proteins.faa"
millard_metadata = "" # path to the Millard phage metadata file, e.g. "/mnt/shared/scratch/vhoikkal/private/databases/millard_phages/4May2024_data.tsv"

temperature_data = "data/200617_TEMPURA.csv" #path to the temperature data file

genomes_folder = "" # path to the folder containing bacterial genomes, e.g. "/home/vhoikkal/scratch/private/databases/ncbi_genomes/bacteria/ncbi_dataset/data"
archaea_folder = "" # path to the folder containing archaeal genomes, e.g. "/home/vhoikkal/scratch/private/databases/ncbi_genomes/archaea/ncbi_dataset/data"
genome_count = 51920 # number of genomes to sample from the genomes folder. This is used only if getGenomesBy is set to "random".
subsampling_seed = 666 # seed for random subsampling of genomes. This is used only if getGenomesBy is set to "random".

genomes_json = os.path.join(genomes_folder,"assembly_data_report.jsonl")

list_of_RNs = ["crn1", "crn2", "crn3", "csx15", "csx16", "csx20", "crn4a", "crn4b"]

if cas10_anchor == False:
    prefiltering_wildcards = "05_host_genomes"
    prefiltering_host_genomes = "06_host_genomes"
elif cas10_anchor == True:
    prefiltering_wildcards = "02_host_wildcards"
    prefiltering_host_genomes = "03_host_genomes"

def aggregate_crisprcas(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/07_cctyper/{i}/{i}_renaming.done", i=ivals)

def aggregate_host_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/02_genome_wildcards/{i}.txt", i=ivals)

def aggregate_cas10_sequences(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa", c=cvals)


def aggregate_renamed_crisprs(wildcards):
    checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
    ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    return expand(base_path + "/07_cctyper/{j}/{j}_renaming.done", j=ivals)

def aggregate_download_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{i}/{i}_genome.fna", i=ivals)


    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    print(expand(base_path + "/19_16S_seq/{j}_16S.fna", j=ivals))
    return expand(base_path + "/19_16S_seq/{j}_16S.fna", j=ivals)

def aggregate_cas10_booleans(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/062_genomes_cas10/{i}/{i}_cas10_boolean.tsv", i=ivals)

def aggregate_cas10_seq_postboolean(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/063_cas10_seq/{i}/{i}_cas10.faa", i=ivals)


def aggregate_cas10_sequences_prior_to_1st_clustering(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/063_cas10_seq/{i}/{i}_cas10.faa", i=ivals)


def aggregate_cas10_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/03_postfiltering_genome_wildcards/{j}.txt", j=jvals)


    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/30_clustertable/{j}/{j}.tsv", j=jvals)

def aggregate_typeIII_info(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_crispr_iii_info.tsv", c=cvals)

def aggregate_taxInfo(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{j}/{j}_taxon.txt", j=ivals)

def aggregate_renamed(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/07_cctyper/{j}/{j}_renaming.done", j=ivals)



    checkpoint_output = checkpoints.align_matcher.get(**wildcards).output[0]
    kvals = glob_wildcards(os.path.join(checkpoint_output,"{k}.done")).k
    return expand(base_path + "/42_matching_trees_distance_matrices/{k}/{k}_correlation.tsv", k=kvals)

def aggregate_cATyper_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv", c=cvals)

def aggregate_ring_nuclease_fusions(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/73_ring_nuclease_fusions/{c}/{c}_ring_nuclease_fusions.tsv", c=cvals)


def aggregate_validated_new_effectors_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv", c=cvals)

def aggregate_validated_new_effectors_analysis(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_ring_nucleases_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv", c=cvals)

def aggregate_ring_nucleases_analysis(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_ring_nucleases_analysis_no_length_lim(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_ring_nucleases_cas10_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/101_ring_fusions_cas10/{c}/{c}_ring_nuclease_cas10_hmm.tsv", c=cvals)

def aggregate_cATyper_analysis(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_cATyper_etp(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_effector_to_protein.tsv", c=cvals)

def aggregate_cATyper_pte(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_protein_to_effector.tsv", c=cvals)

def aggregate_validated_new_effectors_etp(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_effector_to_protein.tsv", c=cvals)

def aggregate_validated_new_effectors_pte(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_protein_to_effector.tsv", c=cvals)

def aggregate_ring_nucleases_etp(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis/{c}/{c}_effector_to_protein.tsv", c=cvals)

def aggregate_ring_nucleases_pte(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis/{c}/{c}_protein_to_effector.tsv", c=cvals)

def aggregate_ring_nucleases_pte_no_length_lim(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}/{c}_protein_to_effector.tsv", c=cvals)

def aggregate_ring_nucleases_etp_no_length_lim(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}/{c}_effector_to_protein.tsv", c=cvals)


def aggregate_CorA_sequences_catyper(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    #if a file does not exist, touch it
    list_of_coras = []
    for c in cvals:
        #check if this locus has a cas10.faa file. We don't want to study CorAs that are not associated with Cas10less loci or those
        #in which filtering removed Cas10 due to length
        #check if cas10 file exists
        #print("Size of cas10 file for locus " + c + " is " + str(os.path.getsize(base_path + "/09_crispr_iii_CorA/loci/" + c + "/" + c + "_Cas10.faa"))) 
        #if directory exists
        cas10_dir_path = base_path + "/09_crispr_iii_CorA"
        cora_dir_path = base_path + "/11_cATyper_analysis/" + c
        if (os.path.exists(cas10_dir_path)) & (os.path.exists(cora_dir_path)):
            if os.path.getsize(base_path + "/09_crispr_iii_CorA/loci/" + c + "/" + c + "_Cas10.faa") > 0:
                #print("Cas10 file is large enough for locus " + c + " so we will add it to the list of coras")
                if c == "GCF_003544875.1_0":
                    print("This is the weird one")
                list_of_coras.append(c)
                if not os.path.exists(base_path + "/11_cATyper_analysis/" + c + "/cora.faa"):
                    print("Touching an empty CorA")
                    open(base_path + "/11_cATyper_analysis/" + c + "/cora.faa", 'a').close()
            if c == "GCF_003544875.1_0":
                print("This is the weird one 2")
    print("Length of list of coras: " + str(len(list_of_coras)))
    #return the list of coras after adding the filepaths to them
    fullpath_coras = []
    for locus in list_of_coras:
        fullpath_coras.append(base_path + "/11_cATyper_analysis/" + locus + "/cora.faa")
    print(fullpath_coras)
    return fullpath_coras

def aggregate_CorA_sequences(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/CorA.faa", c=cvals)


def aggregate_unknowns(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/30_unknown_effectors/{c}/{c}_unknown_proteins.faa", c=cvals)

def aggregate_unknowns_locus_info(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/30_unknown_effectors/{c}/{c}_locus_unknown_info.tsv", c=cvals)

def aggregate_crispr_locus_proteins(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/072_crispr_locus_proteins/{c}/{c}_crispr_locus_proteins.faa", c=cvals)

def aggregate_locus_viz(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/90_locus_viz/{c}/{c}_viz.png", c=cvals)

def aggregate_group4_blasts(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/group4_prober/{c}/{c}.blast", c=cvals)

def aggregate_known_effector_wildcards(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/45_effector_tree/{effector}_tree.txt", effector=effector_vals)

def aggregate_known_effector_wildcarder(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/40_known_effector_wildcards/{effector}.eff", effector=effector_vals)

def aggregate_known_effector_annotations(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered_mapped.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite/{effector}/{effector}_all.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite_parser(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite/{effector}/{effector}_hhsuite_parsed.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite_parser_cogs(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_hhsuite_parsed_cogs.tsv", effector=effector_vals)

def aggregate_cora_neighbourhoods(wildcards):
    checkpoint_output = checkpoints.cora_neighbourhood_preparation.get(**wildcards).output[0]
    cora_loci = glob_wildcards(os.path.join(checkpoint_output,"{cora_locus}_crispr_locus_proteins.faa")).cora_locus
    return expand(base_path + "/52_cora_neighbour_analysis/{cora_locus}/neighbourhood_results.tsv", cora_locus=cora_loci)

#Cas10 blast parameters
blast_headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sacc"
blast_headers_group4 = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tlocus"
cas10_e_value = "1e-20"
corA_hmm_value = "1e-20"

group4_pfam_headers = "target_name\taccession\tquery_name\taccession\tE-value\tscore\tbias\tE-value\tscore\tbias\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tdescription_of_target"

print("Starting ring nuclease / transcription factor pipeline (Crn4 edition)")

rule all: 
    input: base_path + "/done"
    #input: base_path + "/01_genomelist/annotated_genomelist.csv"
    #input: base_path + "/062_genomes_cas10/{i}/{i}_cas10.tsv"
    #input: base_path + "/064_cas10_clusters/cas10_all.faa"
     


rule write_down_genome_info:
    '''
    Outputs a text file genome names based on folders in genomes_folders
    '''
    output: base_path + "/01_genomelist/genomelist.txt"
    threads: thread_ultrasmall
    shell:
        """
        cd {genomes_folder}
        find * -maxdepth 0 -type d > {output}
        """

rule annotate_bacteria_and_archaea_domains:
    """
    Marks down whether a sample is an archaeon or a bacterium
    1. Create a .csv file in the **archaea folder**. This two-column file shows that each sample ID there is, is an archaeon
	2. Create another .csv file in the **bacteria folder**. Note that bacteria folder also contains the archaea (they have been copied there to simplify the fetching of genomes).
        This csv will claim that all samples, including the archaea, are bacteria.
	3. Use the archaea-csv file to **reannotate** the .csv file created in step two, making all archaea appear archaea
	4. The final output file will then be the one created in step three

    This table can then be used in subsequent rules by merging it to any output table using sample name as common identifier.
    *Note that this solution is not ideal and will break down easily if the input genomes are in a different format etc.
    So make sure that the bacteria folder contains **all** samples and that the archaea folder only contains archaeons*
    """
    output: base_path + "/01_genomelist/annotated_genomelist.csv"
    params: 
        archaea_csv = base_path + "/01_genomelist/archaea.csv",
        bacteria_csv = base_path + "/01_genomelist/bacteria.csv"
    threads: thread_ultrasmall
    shell:
        """
        cd {archaea_folder}
        find * -maxdepth 0 -type d > {params.archaea_csv}

        cd {genomes_folder}
        find * -maxdepth 0 -type d > {params.bacteria_csv}

        cd {program_root}
        python3 scripts/annotate_bacteria_and_archaea.py --archaea {params.archaea_csv} --bacteria {params.bacteria_csv} --out {output}
        """


if getGenomesBy == "local":
    checkpoint expand_host_genomelist_to_wildcards:
        '''
        Takes as input a list of genome accessions separated by newline. These genomes are then turned into wildcards.
        In random mode, uses subsampling to randomly sample a set of genomes.
        '''
        input:
            genome_info = rules.write_down_genome_info.output,
            domain_annotations = rules.annotate_bacteria_and_archaea_domains.output
        output: directory(base_path + "/" + prefiltering_wildcards)
        threads: thread_ultrasmall
        run:
            import pandas as pd
            if not os.path.exists(str(output)):
                os.makedirs(str(output))
            genomefiles = pd.read_csv(str(input.genome_info), sep = "\t", header = None)
            #print("HELLOO " + genomefiles)
            genomes = [] #final list of accessions

            #Random mode. Subsamples n genomes from all genomes randomly with seed.
            if config["genome_mode"] == "random":
                subsample_genomes = genomefiles.sample(n = genome_count, random_state=subsampling_seed)
                genomes = subsample_genomes[0]

            elif config["genome_mode"] == "all":
                genomes = genomefiles[0]

            #Taxid id mode. NCBI datasets json is probed to find representatives of a taxid.
            elif config["genome_mode"] == "taxid":
                json_df = pd.read_json(genomes_json, lines=True)
                json_df = json_df.join(pd.json_normalize(json_df['organism'])).drop('organism', axis='columns') #expand pandas dictionary column into new columns in the actual df
                taxidlistlist = []

                with open(config["taxidlistfile"]) as f: #read taxid list file
                    taxids_comma = f.read()
                    idlist = taxids_comma.split(",")

                taxidlist = pd.DataFrame(idlist) #convert list to df
                #taxidlist = pd.read_csv(config["taxidlistfile"], sep=",", header=None)
                taxidlist.columns = ["taxId"] #add column to the df
                json_df['taxId'] = json_df['taxId'].astype("string")
                taxidlist['taxId'] = taxidlist['taxId'].astype("string")
                chosen = pd.merge(json_df, taxidlist, on="taxId", how="inner")
                #chosen = json_df.loc[json_df["taxId"] == int(config["taxid"])]
                print("Shape of subsampled dataframe using TaxIDs from " + str(config["taxidlistfile"]) + ": " + str(chosen.shape))
                genomes = chosen["accession"].values.tolist()


            #add genome to list if .faa, .gff and .fna files are present
            for i in genomes:
                if (os.path.isfile(os.path.join(genomes_folder,i,"protein.faa")) and (os.path.isfile(os.path.join(genomes_folder,i,"genomic.gff")))): #make sure it's annotated
                    for fname in os.listdir(os.path.join(genomes_folder,i)): #checking fna is more complicated because the filename is not consistent
                        if fname.endswith('.fna'):
                            with open(str(output) + "/" + str(i) + '.txt', 'w') as f:
                                f.write(str(output))
        


    rule download_genomes:
        '''
        Based on the wildcards from expand_host_genomelist_to_wildcards,
        makes symlinks for the desired genomes from the genomes_folder
        into the working folder (prefiltering_host_genomes).
        '''
        input:
            ivalue = base_path + "/" + prefiltering_wildcards + "/{i}.txt",
        output: 
            fna = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_genome.fna",
            faa = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa",
            gff = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_features.gff"
        params: #these are mostly deprecated
            folder = base_path + "/" + prefiltering_host_genomes,
        conda: "envs/ncbidownload.yaml"
        threads: thread_ultrasmall
        log:
            out = base_path + "/logs/" + prefiltering_host_genomes + "/{i}.log",
            err = base_path + "/logs/" + prefiltering_host_genomes + "/{i}.err"
        shell:
            '''
            echo "Creating symlinks for chosen genomes"
            ln -s {genomes_folder}/{wildcards.i}/*.fna {output.fna}
            ln -s {genomes_folder}/{wildcards.i}/*.faa {output.faa}
            ln -s {genomes_folder}/{wildcards.i}/*.gff {output.gff}
            '''


elif getGenomesBy == "remote":
    species_space = str(config["species"])
    checkpoint getHosts:
        '''
        Downloads desired host genomes from NCBI using the wonderful python program ncbi-genome-download.
        Files are then gunzipped and named after their parent folder
        '''
        output: directory(base_path + "/" + prefiltering_host_genomes)
        conda: "envs/ncbidownload.yaml"
        params:
            parallels = thread_hogger,
            folder = base_path + "/06_host_genomes"
        threads: thread_ultrasmall
        shell:
            '''
            mkdir -p {params.folder}
            cd {params.folder}
            count=$(ncbi-genome-download --dry-run --genera "{species_space}" bacteria | wc -l)
            echo "Will download $count genomes from {species_space}. Note that the download will likely fail if the dataset is large (>1000 genomes). In such a case, just restart the script."
            printf 'Do you want to continue? (y/n)? '
            read -p "Do you want to continue? (y/n) " -n 1 -r
            if [[ $REPLY =~ ^[Yy]$ ]] ;then
                ncbi-genome-download --genera "{species_space}" bacteria --parallel {params.parallels} --formats fasta,protein-fasta,gff,assembly-report
                mv refseq/bacteria/* {params.folder}
                find {params.folder} -type f -name '*.gz' | tqdm | xargs gunzip
                echo renaming
                rename 's!(.*)/(.*)genomic\.fna!$1/$1_genome.fna!' */*genomic.fna
                rename 's!(.*)/(.*)protein\.faa!$1/$1_proteins.faa!' */*protein.faa
                rename 's!(.*)/(.*)genomic\.gff!$1/$1_features.gff!' */*genomic.gff
                rename 's!(.*)/(.*)assembly_report\.txt!$1/$1_report.txt!' */*assembly_report.txt
            else
                exit 1
            fi

            '''


if cas10_anchor == True: #if we filter analyzable genomes by the presence of Cas10. Uses Cas10 HMM profiles from CCtyper.
    print("Running Cas10 stuff")
    rule Cas10_genomes:
        '''
        First filter all sequences to get only >500 AA long proteins. This prevents the inclusion of
        truncated Cas10s or those that are cut artifically due to contig ending. Also reduces HMM searches.
        
        Then:
        Searches the filtered proteins of host against a user-specified preprepared hmm profile.
        Strips file of "#" -rows and checks if there is anything left. If there is, the gene has been found.
        Boolean output format:
        i,True/False,cas10_accession
        '''
        input:
            proteins = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa"
        output:
            hmm = base_path + "/062_genomes_cas10/{i}/{i}_cas10.tsv",
            boolean = base_path + "/062_genomes_cas10/{i}/{i}_cas10_boolean.tsv"
        params:
            proteins_filt = base_path + "/062_genomes_cas10/{i}/{i}_proteins_lengthfiltered.faa",
            out = base_path + "/062_genomes_cas10/{i}/{i}_temp.out",
            rows1 = base_path + "/062_genomes_cas10/{i}/{i}_temp_rows_1.out",
            cas10_db = cas10_db,
            rows = base_path + "/062_genomes_cas10/{i}/{i}_temp_rows.out",
            headers = base_path + "/062_genomes_cas10/{i}/{i}_headers.out",
            all_data = base_path + "/062_genomes_cas10/{i}/{i}_all_data.out",
        conda: "envs/hmmer.yaml"
        threads: thread_ultrasmall
        log:
            out = base_path + "/logs/062_genomes_cas10/{i}.log",
            err = base_path + "/logs/062_genomes_cas10/{i}.err"
        shell:
            '''
            echo "Running rule Cas10_genomes"
            echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm}
            if [ -s "{input.proteins}" ]; then
                cat {input.proteins} | seqkit seq -m 500 > {params.proteins_filt}
                if [ -s "{params.proteins_filt}" ]; then
                    hmmscan --domtblout {params.out} --cpu {threads} -E {corA_hmm_value} {params.cas10_db} {params.proteins_filt} 2> {log.err} 1> {log.out}
                    grep -v "#" {params.out} > {params.rows}||:
                    head -1 {params.rows} > {params.rows1}
                    echo "id,cas10_boolean,cas10_acc" > {output.boolean}
                    if [ -s {params.rows1} ]; then
                        cat {params.rows1} >> {output.hmm}
                        ACC=$(awk -F ' ' '{{print $4}}' {params.rows1})
                        echo "{wildcards.i},True","${{ACC}}" > {output.boolean}
                    else
                        echo "{wildcards.i},False,-" > {output.boolean}
                    fi
                    touch {output.hmm}
                else
                    echo "{wildcards.i},False,-" > {output.boolean}
                fi
            else
                echo "{wildcards.i},False,-" > {output.boolean}
            fi
            '''

        
    rule extract_cas10_sequence:
        input:
            boolean = rules.Cas10_genomes.output.boolean,
            proteins = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa"
        output: base_path + "/063_cas10_seq/{i}/{i}_cas10.faa"
        threads: thread_ultrasmall
        conda: "envs/hmmer.yaml"
        shell: 
            '''
            CAS10_ID=$(awk -F ',' '{{print $3}}' {input.boolean})
            if [ "$CAS10_ID" != "-" ]; then
                echo ">"{wildcards.i} > {output}
                seqkit grep -r -p ${{CAS10_ID}} {input.proteins} | seqkit seq -w 0 | tail -n +2 >> {output}
            else
                touch {output}
            fi
            '''

    rule concat_cluster_cas10_proteins:
        '''
        Takes all Cas10 sequences, concatenates them into one fasta and
        clusters this fasta at cas10_cluster_threshold (0-1) similarity.
        Greps only representative lines "marked by *"
        A python script extracts the accession from these lines into output.reps
        '''
        input: aggregate_cas10_seq_postboolean
        output:
            all_cas10 = base_path + "/064_cas10_clusters/cas10_all.faa",
            clusters = base_path + "/064_cas10_clusters/cas10_clust.faa.clstr",
            proteins = base_path + "/064_cas10_clusters/cas10_clust.faa",
            reps = base_path + "/064_cas10_clusters/cas10_unique_genomes.txt",
        params:
            clusterlines = base_path + "/064_cas10_clusters/clusterlines.txt"
        threads: thread_hogger
        conda: "envs/groupChar.yaml"
        shell:
            '''
            echo "concatenating cas10 sequences for clustering"
            find '{base_path}/063_cas10_seq' -maxdepth 2 -type f -wholename '*/*_cas10.faa' -print0 | xargs -0 cat > {output.all_cas10}
            echo clustering
            cd-hit -i {output.all_cas10} -o {output.proteins} -c {cas10_cluster_threshold} -n 5 -d 0 -M 16000 -T {threads}
            grep "*" {output.clusters} > {params.clusterlines}
            python scripts/getClusterRep.py --input {params.clusterlines} --output {output.reps}
            '''

    checkpoint expand_host_genomelist_to_wildcards_postcas10:
        '''
        '''
        input:
            genomes = rules.concat_cluster_cas10_proteins.output.reps
        output: directory(base_path + "/03_postfiltering_genome_wildcards")
        threads: thread_ultrasmall
        run:
            print("Expanding host list to wildcards...")

            list_of_genomes_to_remove = ["GCA_019977735.1"] #use this list to manually exclude any unwanted genomes

            if not os.path.exists(str(output)):
                os.makedirs(str(output))

            with open(str(input.genomes)) as filehandle:
                genomes = [line.rstrip('\n') for line in filehandle]

            for i in genomes:
                if i not in list_of_genomes_to_remove:
                    sample = str(str(i).split(",")[0])
                    with open(str(output) + "/" + sample + '.txt', 'w') as f:
                        f.write(sample)



    rule download_genomes_postCas10:
        '''
        After filtering by presence of Cas10, create new symlinks for such genomes.
        '''
        input: 
            genomes = base_path + "/03_postfiltering_genome_wildcards/{j}.txt",
        output: 
            fna = base_path + "/06_host_genomes/{j}/{j}_genome.fna",
            faa = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
            gff = base_path + "/06_host_genomes/{j}/{j}_features.gff"
        params: #these are mostly deprecated
            folder = base_path + "/06_host_genomes",
            zipped = base_path + "/06_host_genomes/{j}.zip",
            original_fna = base_path + "/06_host_genomes/{j}/*.fna", #for checking successful download
            original_prot = base_path + "/06_host_genomes/{j}/protein.faa",#for checking successful download
            original_gff = base_path + "/06_host_genomes/{j}/genomic.gff", #for checking successful download
        conda: "envs/ncbidownload.yaml"
        threads: thread_ultrasmall
        log:
            out = base_path + "/logs/06_host_genomes/{j}.log",
            err = base_path + "/logs/06_host_genomes/{j}.err"
        shell:
            '''
            ln -s {genomes_folder}/{wildcards.j}/*.fna {output.fna}
            ln -s {genomes_folder}/{wildcards.j}/*.faa {output.faa}
            ln -s {genomes_folder}/{wildcards.j}/*.gff {output.gff}
            '''

rule getTaxInfo:
    input: base_path + "/03_postfiltering_genome_wildcards/{j}.txt"
    output:
        taxon = base_path + "/06_host_genomes/{j}/{j}_taxon.txt"
    threads: thread_ultrasmall
    run:
        import pandas as pd
        json_df = pd.read_json(genomes_json, lines=True)
        json_df = json_df.join(pd.json_normalize(json_df['organism'])).drop('organism', axis='columns') #expand pandas dictionary column into new columns in the actual df
        json_df = json_df.join(pd.json_normalize(json_df['averageNucleotideIdentity'])).drop('averageNucleotideIdentity', axis='columns') #expand pandas dictionary column into new columns in the actual df 
        try:
            row = json_df.loc[json_df["accession"] == wildcards.j]
            species = row['submittedSpecies'].iloc[0].strip("][")
            genus = species.split(" ")[0]
            with open(output.taxon, "w") as f:
                f.write("Sample\tGenus\tSpecies\n")
                f.write(wildcards.j + "\t" + genus + "\t" + species + "\n")
        except:
            with open(output.taxon, "w") as f:
                            f.write(wildcards.j + "\tUnknown\tUnknown\n")
            pass


rule concat_taxInfo:
    input: aggregate_taxInfo
    output: base_path + "/06_host_genomes/taxInfo.txt"
    threads: thread_small
    shell:
        """
        echo "Sample\tgenus\tspecies\n" > {output}
        find '{base_path}/06_host_genomes/' -maxdepth 2 -type f -wholename '*_taxon.txt' -print0 | xargs -0 cat >> {output}
        """

rule CRISPRCasTyper:
    '''
    Types CRISPR-Cas loci in the genomes based on fasta.
    Takes output of the checkpoint. Automatically scales to the number of files
    produced by the checkpoint, so no need to call for the aggregator function here.

    '''
    input: base_path + "/06_host_genomes/{j}/{j}_genome.fna"
    output: base_path + "/07_cctyper/{j}/{j}.done"
    conda: "envs/CRISPRCasTyper.yaml"
    log: base_path + "/logs/07_cctyper/{j}/{j}_cctyper.log"
    params:
        outdir = base_path + "/07_cctyper/{j}"
    threads: thread_ultrasmall
    shell:
        '''
        rm -rf {params.outdir}
        cctyper '{input}' '{base_path}/07_cctyper/{wildcards.j}' --prodigal single 2>{log}
        touch '{output}'
        '''

rule CRISPRCasTyper_rename:
    '''
    UPDATED VERSION. This uses the CRISPR-Cas.tab file to include only fully functional CRISPR-Cas loci
    The sed adds assembly information to the results file, if it exists (which is originally based on just the nucleotide accessions)
    '''
    input: rules.CRISPRCasTyper.output
    output: base_path + "/07_cctyper/{j}/{j}_renaming.done"
    params:
        outdir = base_path + "/07_cctyper/{j}"
    threads: thread_ultrasmall
    shell:
        '''
        if test -e "{params.outdir}/CRISPR_Cas.tab";then
            mv {params.outdir}/CRISPR_Cas.tab {params.outdir}/CRISPR_Cas_old.tab
            sed '1s/$/\tassembly/; 2,$s/$/\t{wildcards.j}/' {params.outdir}/CRISPR_Cas_old.tab > {params.outdir}/CRISPR_Cas.tab
            rm {params.outdir}/CRISPR_Cas_old.tab
        fi
        touch {output}
        '''

rule concat_renamed_crisprs:
    '''
    This rule aggregates the outputs of the CRISPRCasTyper_rename rule.
    It is necessary to do this because the CRISPRCasTyper_rename rule is called by the checkpoint,
    and the checkpoint does not automatically aggregate the outputs of the rule.
    '''
    input: aggregate_renamed
    output: base_path + "/07_cctyper/renaming.done"
    threads: thread_ultrasmall
    shell:
        '''
        touch {output}
        ''' 

checkpoint type_iii_wildcarder:
    '''
    Extracts only type III loci from the CCTyper output
    This checkpoint examines outputs from CCTyper and creates new wildcards based on complete loci.
    This marks a point where the pipeline no longer looks at individual strains, but rather at the CRISPR-Cas loci individually.
    As outputs, this checkpoint extracts locus-specific information from CCTyper outputs (CRISPR-Cas.tab and cas_operons.tab)
    NOTE: now excludes all loci with >1 Cas10s
    '''
    input:
        cctyper_done = aggregate_renamed,
    output: directory(base_path + "/071_cctyper_loci")
    params:
        interference_cutoff = crispr_locus_interference_cutoff, #0-100. If the interference score is lower than this, the locus is discarded
        cctyper_folder = base_path + "/07_cctyper",
    threads: thread_hogger
    shell:
        '''
        python3 scripts/loci_wildcarder.py --input_folder {params.cctyper_folder} --output_folder {output} --interference_cutoff {params.interference_cutoff}
        '''


rule cATyper_hmm_search:
    '''
    Runs HMMscan for all proteins in each locus. As database we use the concatenated hmm database of all known effectors.
    The main output file is hmm_rows, which contains all hmm hits for all proteins in the locus.
    The column target_name is the hmm hit and query_name is the protein accession
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
    output:
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        temp_rows = base_path + "/10_cATyper_hmm/{c}/{c}_temp_rows.out",
        hmm_rows = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv",
        cora = base_path + "/10_cATyper_hmm/{c}/CorA.faa",
    params:
        outdir = base_path + "/10_cATyper_hmm/{c}",
        hmm_msa_folder = hmm_msa_folder,
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_profile = hmm_database_folder + "/" + hmm_database_file,
        temp_hmm = base_path + "/10_cATyper_hmm/{c}/{c}_temp.out",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/10_cATyper_hmm/logs/{c}.out",
        err = base_path + "/10_cATyper_hmm/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode pre_hmm 2> {log.err} 1> {log.out}
        echo "Running hmmscan" >> {log.out}
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {catyper_hmm_evalue} {params.hmm_profile} {output.contig_proteins}  &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        touch {output.cora}
        '''

rule concatenate_cATyper_hmm:
    input: aggregate_cATyper_hmm
    output: base_path + "/10_cATyper_hmm/cATyper_all.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule cATyper_analysis:
    '''
    This is the analysis step of cATyper. Note that we are still using the same python script as in the previous rule,
    but the mode is now set to "post_hmm".
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv",
    output:
        catyper = base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/11_cATyper_analysis/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/11_cATyper_analysis/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/11_cATyper_analysis/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/11_cATyper_analysis/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = hmm_msa_folder,
        tmhmm_model_path = TM_path + "/TMHMM2.0.model",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/11_cATyper_analysis/logs/{c}.out",
        err = base_path + "/11_cATyper_analysis/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        echo "Running cATyper analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "known_effector" --effector_plot_data {output.plottable_effector_positions} --tmhmm_model_path {params.tmhmm_model_path} 2> {log.err} 1> {log.out}
        '''

rule concatenate_cATyper_analysis:
    input: aggregate_cATyper_analysis
    output: base_path + "/11_cATyper_analysis/cATyper_all.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {base_path}/11_cATyper_analysis/*/*_cATyper_results.tsv > {output}
        touch {output}
        """

rule concatenate_cATyper_analysis_effector_scores:
    '''
    Concatenates info on the effector hits from cATyper analysis
    '''
    input:
        protein_to_effector = aggregate_cATyper_pte,
        effector_to_protein = aggregate_cATyper_etp
    output:
        protein_to_effector_concatenated = base_path + "/11_cATyper_analysis/cATyper_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/11_cATyper_analysis/cATyper_effector_to_protein.tsv",
    threads: thread_small
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule analyse_cATyper_effector_scores:
    '''
    The algorithm that analyses how well each hmm profile performs in cATyper.
    '''
    input:
        pte = rules.concatenate_cATyper_analysis_effector_scores.output.protein_to_effector_concatenated,
        etp = rules.concatenate_cATyper_analysis_effector_scores.output.effector_to_protein_concatenated,
    output:
        effector_scores = base_path + "/11_cATyper_analysis/cATyper_effector_scores.tsv",
        effector_scores_summary = base_path + "/11_cATyper_analysis/cATyper_effector_scores_summary.tsv",
        #effector_scores_plot1 = base_path + "/11_cATyper_analysis/cATyper_effector_scores_plot1.png",
    params:
        outdir = base_path + "/11_cATyper_analysis"
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    shell:
        '''
        python scripts/catyper_effector_scores.py --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        touch {output.effector_scores}
        touch {output.effector_scores_summary}
        '''

checkpoint effector_wildcarder:
    '''
    This rule creates wildcards for each previosly known effector.
    Wildcards are used downstream when creating phylogenetic trees etc for each effector.
    '''
    output:
        directory(base_path + "/40_known_effector_wildcards")
        #hmm_targets = base_path + "/40_known_effector_wildcards/{effector}/{effector}.txt"
    params:
        hmm_targets = base_path + "/40_known_effector_wildcards/alltargets.txt",
        hmm_msa_folder = hmm_msa_folder,
        outdir = base_path + "/40_known_effector_wildcards",
    threads: thread_small
    run:
        import glob
        #Equivalent of `mkdir`
        Path(params.outdir).mkdir(parents=True, exist_ok=True)
        
        #Equivalent of `ls {hmm_msa_folder}/*/*.hmm > {hmm_targets}`
        with open(params.hmm_targets, 'w') as f:
            files = glob.glob(f"{params.hmm_msa_folder}/*/*.hmm") #get all files in the hmm_msa_folder
            f.write('\n'.join(files))

        with open(params.hmm_targets, 'r') as f:
            lines = f.readlines()
            for line in lines:
                filename = os.path.basename(line) # get filename from path
                effector = re.split("\.", filename)[0] #remove extension
                if "#" not in line:
                    effector = re.split("_", effector)[1]
                elif "#" in line: #in case a single effector is characterised by multiple hmm profiles, we mark the filenames with # followed by the effector, e.g ca6_csm6-ca6_italicus#csm6-ca6
                    effector = re.split("#", effector)[1]
                print(effector)
                with open(f'{params.outdir}/{effector}.eff', 'w') as f:
                    f.write(effector)

rule effector_fetcher:
    '''
    Fetches protein sequences for an effector from catyper_analysis outputs
    '''
    input:
        aggregate_known_effector_wildcarder,
        catyper_done = rules.concatenate_cATyper_analysis.output, #this indicates we can safely probe the folders for effector sequences
    output:
        multifasta = base_path + "/41_known_effector_mf/{effector}.faa",
    threads: thread_ultrasmall
    shell:
        '''
        effector_name={wildcards.effector}
        echo "Fetching effector sequences for {wildcards.effector}"
        touch {output.multifasta}
        find {base_path}/11_cATyper_analysis/ -name \"*${{effector_name}}*.faa\" -exec cat {{}} \; >> {output.multifasta}
        '''

rule effector_hmmer:
    '''
    Uses hmmscan against pfam to annotate the discovered effectors' domains.
    Output used when plotting tree with R
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        raw_hmm = base_path + "/42_effector_hmmer/{effector}.hmm",
        #filtered_hmm_tsv = base_path + "/42_effector_hmmer/{effector}.tsv",
    params:
        base_path = base_path + "/42_effector_hmmer",
        pfam_db = pfam_db
    threads: thread_small
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/42_effector_hmmer/logs/{effector}.out",
        err = base_path + "/42_effector_hmmer/logs/{effector}.err",
    shell:
        '''
        python3 scripts/effector_hmmer.py --input {input.multifasta} --output_basepath {params.base_path} --effector {wildcards.effector} --pfam_path {params.pfam_db} 2> {log.err} 1> {log.out}
        '''

rule cATyper_hmm_search_hhsuite_pdb:
    '''
    HMM search against the PDB database using HHsuite.
    Using a custom made wrapper to first divide the multifasta into separate fastas for each protein as required by hhblits.
    The wrapper divides a multifasta (wildcard {c}) into its constituent proteins and then runs hhblits on each protein.
    The resulting files are .tsv files, each corresponding to a single protein within the effector wildcard (e.g. RelE protein X)
    These individual .tsv files are then concatenated into a single .tsv file for each effector and fed to the parse_hhsuite rule.
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        hhsuite_concat = base_path + "/42_effector_hhsuite/{effector}/{effector}_all.tsv",
    params:
        outdir = base_path + "/42_effector_hhsuite/{effector}",
        pdb30 = pdb30_db + "/pdb30",
    conda: "envs/hhsuite.yaml"
    threads: thread_small
    log:
        out = base_path + "/42_effector_hhsuite/logs/{effector}.out",
        err = base_path + "/42_effector_hhsuite/logs/{effector}.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.pdb30}
        cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        '''

rule parse_hhsuite:
    '''
    Parser output from HHSuite (PDB)
    '''
    input: rules.cATyper_hmm_search_hhsuite_pdb.output.hhsuite_concat
    output: base_path + "/42_effector_hhsuite/{effector}/{effector}_hhsuite_parsed.tsv"
    conda: "envs/hhsuite.yaml"
    threads: thread_small
    params:
        database = "PDB"
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output}  --database {params.database}
        '''

rule concatenate_cATyper_hmm_hhsuite:
    input: aggregate_cATyper_hhsuite_parser
    output: base_path + "/42_effector_hhsuite/cATyper_all.tsv"
    threads: thread_small
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule cATyper_hmm_search_hhsuite_cog:
    '''
    HMM search against the COG database using HHsuite.
    Using a custom made wrapper to first divide the multifasta into separate fastas for each protein as required by hhblits.
    The wrapper divides a multifasta (wildcard {c}) into its constituent proteins and then runs hhblits on each protein.
    The resulting files are .tsv files, each corresponding to a single protein within the effector wildcard (e.g. RelE protein X)
    These individual .tsv files are then concatenated into a single .tsv file for each effector and fed to the parse_hhsuite rule.
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        hhsuite_concat = base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_all.tsv",
    params:
        outdir = base_path + "/42_effector_hhsuite_cogs/{effector}",
        cogs = cogs_db,
    conda: "envs/hhsuite.yaml"
    threads: thread_small
    log:
        out = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.out",
        err = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.cogs} 2> {log.err} 1> {log.out}
        #check if any .tsv files exist in {params.outdir}/hhblits
        if [ -s {params.outdir}/hhblits/*.tsv ]; then
            cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        else
            touch {output.hhsuite_concat}
        fi
        '''

rule parse_hhsuite_cogs:
    input: rules.cATyper_hmm_search_hhsuite_cog.output.hhsuite_concat
    output: base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_hhsuite_parsed_cogs.tsv"
    conda: "envs/hhsuite.yaml"
    params:
        database = "COGs",
        mapping = cogs_db + "/cog-20.def.tab"
    threads: thread_small
    log:
        out = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.out",
        err = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.err",
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output} --database {params.database} --mapping {params.mapping} 2> {log.err} 1> {log.out}
        '''

rule concatenate_cATyper_hmm_hhsuite_cogs:
    input: aggregate_cATyper_hhsuite_parser_cogs
    output: base_path + "/42_effector_hhsuite_cogs/cATyper_all_hmm_cogs.tsv"
    threads: thread_small
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule effector_analyse_domains:
    '''
    Takes the HMM output from effector_hmmer and creates a simplified table of domains for each effector.
    This is only for the Pfam analysis.
    LEGACY: Links protein ID to each effector (the original effector fastas previously contained locus IDs, not protein IDs)
    '''
    input:
        raw_hmm = rules.effector_hmmer.output.raw_hmm,
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        domains = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered_mapped.tsv",
        filtered_hmm_tsv = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered.tsv",
        sorted_hmm_tsv = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted.tsv"
    params:
        base_path = base_path + "/43_effector_hmmer_analysis",
        effector_locus_map = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}.locus_map.tsv",
    threads: thread_small
    log:
        out = base_path + "/43_effector_hmmer_analysis/logs/{effector}.out",
        err = base_path + "/43_effector_hmmer_analysis/logs/{effector}.err"
    shell:
        '''
        cat {base_path}/11_cATyper_analysis/*/*_protein_to_effector.tsv | grep "{wildcards.effector}" > {params.effector_locus_map}
        python3 scripts/effector_hmmer_analyse.py --input {input.raw_hmm} --multifasta {input.multifasta} --output_basepath {params.base_path} --effector {wildcards.effector} --effector_locus_map {params.effector_locus_map} 2> {log.err} 1> {log.out}
        '''

rule effector_align:
    input: rules.effector_fetcher.output.multifasta
    output: base_path + "/44_effector_alignments/{effector}.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/44_effector_alignments/logs/{effector}.out",
        err = base_path + "/44_effector_alignments/logs/{effector}.err"
    threads: thread_small
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out} 
        '''


rule effector_tree:
    input: rules.effector_align.output
    output: base_path + "/45_effector_tree/{effector}_tree.txt",
    conda: "envs/trees.yaml"
    threads: thread_small
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule concatenate_effector_wildcards:
    input:
        trees = aggregate_known_effector_wildcards,
        annotations = aggregate_known_effector_annotations,
        checkpoint = aggregate_known_effector_wildcarder
    output: base_path + "/known_effectors_finished.txt"
    threads: thread_small
    shell:
        '''
        cat {input} > {output}
        '''

rule crispr_locus_proteins:
    '''
    Extracts all proteins from a CRISPR-Cas locus +- 4000 bp
    '''
    input:
        cas_operons_tsv = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
    output:
        crispr_locus_proteins = base_path + "/072_crispr_locus_proteins/{c}/{c}_crispr_locus_proteins.faa",
    params:
        output_folder = base_path + "/072_crispr_locus_proteins/{c}",
        locus = "{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
    threads: thread_ultrasmall
    run:
        import os
        import pandas as pd
        from Bio import SeqIO
        import gffutils
        import re
        import sys

        #when looking for effectors in the locus, the effector search range is the number of bases up or downstream of the cctyper defined cas operon boundaries
        effector_search_range = 4000

        #get sample name by splitting the locus name at the underscore and taking the first two parts
        sample = params.locus.split("_")[0] + "_" + params.locus.split("_")[1]

        #read in the cas operons file
        cas_operons_df = pd.read_csv(input.cas_operons_tsv, sep ="\t", header = 0)

        #find the gff file for the sample from the host_genomes_folder/samples folder
        gff_file = os.path.join(params.host_genomes_folder, sample, sample + "_features.gff")

        #Using gffutils, create a database of the gff file
        db_path = os.path.join(params.output_folder, params.locus)

        #get the path to the proteins fasta file for the current sample
        proteins_fasta = os.path.join(params.host_genomes_folder, sample, sample + "_proteins.faa")

        #create folder for the above path if it does not exist
        #if not os.path.exists(os.path.dirname(db_path)):
        #    print("Creating folder " + str(os.path.dirname(db_path)))
        #    os.makedirs(os.path.dirname(params.db_path))
        db = gffutils.create_db(gff_file, 
                                dbfn=db_path + ".db", 
                                force=True, 
                                keep_order=True, 
                                merge_strategy='merge', 
                                sort_attribute_values=True,
                                id_spec=["protein_id", "Name", "ID"])

        #read in the db
        db = gffutils.FeatureDB(db_path + ".db", keep_order=True)

        #get contig from row 0, column Contig
        contig = cas_operons_df.iloc[0]["Contig"]

        cas_operon_start = int(cas_operons_df['Start'][0]) - effector_search_range
        cas_operon_end = int(cas_operons_df['End'][0]) + effector_search_range

        #from the gff file, extract all the features that are on the contig. The featuretype must be "CDS"
        protein_ids = []
        proteins_on_contig = db.region(featuretype='CDS', seqid=contig)

        #then, extract all proteins whose coordinates are between the start and end of the cas operon +- the effector search range
        proteins_on_crispr_locus = db.region(featuretype='CDS', seqid=contig, start = cas_operon_start, end = cas_operon_end)

        #convert the returned generator to a list
        print("Extracting protein IDs from gff file")
        for i in proteins_on_crispr_locus:
            #check if the attributes of the feature has a key called 'protein_id'
            if "protein_id" in i.attributes:
                id_in_gff = str(i.attributes['protein_id'][0])
                #id = str(i.attributes['ID'][0]).split("-")[1]
                protein_ids.append(id_in_gff)

        #using biopython, extract the protein sequences using the list protein_ids from the proteins fasta file
        #the proteins fasta file is in the same folder as the gff file
        protein_seqs = []
        print("Reading protein fasta file from " + str(proteins_fasta))
        for record in SeqIO.parse(proteins_fasta, "fasta"):
            if record.id in protein_ids:
                protein_seqs.append(record)

        #write the protein sequences to a multifasta file
        print("Writing protein sequences to " + str(output.crispr_locus_proteins))
        SeqIO.write(protein_seqs, output.crispr_locus_proteins, "fasta")



rule aggregate_crispr_locus_proteins:
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    input: aggregate_crispr_locus_proteins
    output: base_path + "/072_crispr_locus_proteins/crispr_locus_proteins_all.faa"
    threads: thread_small
    shell:
        '''
        cat {input} > {output}
        '''



rule typeIII_characterizer:
    '''
    Characterizes previously known genes from type III loci (currently Cas10, Cas7, Cas5, CorA).
    This version incorporates the HD domain search. HD domains are searched for in the Cas10 protein
    and within the first 5 to 30 residues. Also adds Cas10 length.
    Now also looks for GGDD domain in Cas10.
    '''
    input:
        cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv", #the cas_operons file is trimmed down to only include this specific {c} locus
        #gff = base_path + "/06_host_genomes/{j}/{j}_features.gff",
        #proteins = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
    output:
        Cas10_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa",
        Cas5_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas5.faa",
        Cas7_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas7.faa",
        #CorA_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_CorA.faa", deprecated as of 17.1.2024
        info = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_crispr_iii_info.tsv",
    params:
        this_folder = base_path + "/09_crispr_iii_CorA/loci/{c}",
        outputfolder = base_path + "/09_crispr_iii_CorA/loci/{c}",
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        cctyper_folder = base_path + "/07_cctyper"
    conda: "envs/gff_utils.yaml"
    log:
        out = base_path + "/09_crispr_iii_CorA/logs/{c}.out",
        err = base_path + "/09_crispr_iii_CorA/logs/{c}.err"
    threads: thread_ultrasmall
    shell:
        '''
        python3 scripts/type_iii_effector_finder_2.0_HD.py --locus_id {wildcards.c} --sample_folder {params.sample_folder} --this_folder {params.this_folder} --outputfolder {params.outputfolder} --cas_operons {input.cctyper} --info_out {output.info} --cctyper_path {params.cctyper_folder} 2> {log.err} 1> {log.out}
        touch {output.Cas10_fasta}
        touch {output.Cas5_fasta}
        touch {output.Cas7_fasta}
        '''


rule concatenate_type_iii_info:
    '''
    1. Creates header in output
    2. Concatenates individual loci info files
    3. Removes the headers from the concatenated file using inverse grep
    4. Joins this modified file with the header file
    '''
    input: aggregate_typeIII_info
    output: base_path + "/09_crispr_iii_CorA/loci/type_iii_info.tsv"
    params:
        temp_out = base_path + "/09_crispr_iii_CorA/loci/temp.tsv" 
    threads: thread_ultrasmall
    shell:
        '''
        echo "Cas10\tCas5\tCas7\tLocus\tSample\tCas10_GGDD\tCas10_GGDD_coord\tCas10_GGDD_seq\tCas10_GGDE\tCas10_GGDE_coord\tCas10_GGDE_seq\tCas10_GGED\tCas10_GGED_coord\tCas10_GGED_seq\tCas10_HD\tCas10_HD_list\tCas10_DH\tCas10_HD_coord\tCas10_DH_coord\tCas10_coord\tCas10_length\tSubtype\tcyclase_literal" > {params.temp_out}
        find '{base_path}/09_crispr_iii_CorA/loci' -maxdepth 2 -type f -wholename '*/*_crispr_iii_info.tsv' -print0 | xargs -0 tail -n +2 >> {params.temp_out}
        grep -v "==>" {params.temp_out} > {output}
        '''

rule Cas10_concatenate:
    '''
    Concatenates Cas10 sequences output by rule concatenate_type_iii_info
    '''
    input: aggregate_cas10_sequences
    output: base_path + "/10_Cas10_cluster/cas10s.faa"
    threads: thread_ultrasmall
    shell:
        '''
        cat {input} > {output} 
        '''

rule Cas10_HD_hmm_maker:
    '''
    1. Extracts the first 10-35 AA of Cas10s that are marked as having an HD domain.
    2. Aligns the Cas10s using Muscle
    3. Creates a new HD-HMM profile from this alignment

    Also has the option to use a premade alignment. Use --use_existing_alignment to do this (does not need argument)
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        cas10_sequences = rules.Cas10_concatenate.output,
    output:
        hmm = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.hmm",
        msa = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.msa",
        faa = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.faa",
    conda: "envs/hmmer.yaml"
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/Cas10_HD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm} --msa {output.msa} --faa {output.faa}  --existing_alignment_path {modified_cas10_hd}
        touch {output.msa}
        touch {output.faa}
        '''

rule Cas10_HD_hmmer:
    '''
    Using the HD-HMM profiles made in rule Cas10_HD_hmm_maker, searches for HD domains in all Cas10s
    '''
    input:
        hmm_db = rules.Cas10_HD_hmm_maker.output.hmm,
        cas10_sequences = rules.Cas10_concatenate.output
    output:
        raw_table = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits.out"
    params:
        out = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_temp.out",
        rows = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_rows.out",
        rows1 = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_rows1.out",
        E = "1e-01"
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    shell:
        '''
        hmmscan --tblout {params.out} --cpu {threads} -E {params.E} {input.hmm_db} {input.cas10_sequences} &> /dev/null
        grep -v "#" {params.out} > {params.rows}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table}
        cat {params.rows} >> {output.raw_table}
        '''

rule merge_HD_hmm_with_info:
    '''
    Merges the HD-HMM hits with the type III info table
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        HD_hmm_hits = rules.Cas10_HD_hmmer.output.raw_table
    output:
        merged_table = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_merged.tsv"
    threads: thread_ultrasmall
    run:
        import pandas as pd
        info = pd.read_csv(str(input.info_table), sep = "\t")
        print(info)
        hmm_hits = pd.read_csv(str(input.HD_hmm_hits), sep = '\s+')
        print(hmm_hits)

        #reduce the hmm_hits to only query_name and E-value_best_domain
        hmm_hits = hmm_hits[["query_name", "E-value_best_domain"]]

        #rename the columns
        hmm_hits.columns = ["Locus", "HD_E-value"]
        #add column "HD_hmm_boolean" and set it to True
        hmm_hits["HD_hmm_boolean"] = True

        #merge the info table with the hmm_hits table
        merged = pd.merge(info, hmm_hits, on = "Locus", how = "left")
        #if some rows were left without a hit, set the HD_hmm_boolean to False
        merged["HD_hmm_boolean"] = merged["HD_hmm_boolean"].fillna(False)

        #from the final file, remove all columns but Locus, HD_E-value and HD_hmm_boolean
        merged = merged[["Locus", "HD_E-value", "HD_hmm_boolean"]]

        #save the merged table
        merged.to_csv(str(output.merged_table), sep = "\t", index = False)


rule R_HD:
    '''
    R script that produces a histogram of Cas10s and their HD domains
    '''
    input: rules.concatenate_type_iii_info.output
    output:
        HD_histogram = base_path + "/4_R_HD/cas10_HD_lengths.png"
    log:
        out = base_path + "/4_R_HD/logs/HD_hist.out",
        err = base_path + "/4_R_HD/logs/HD_hist.err"
    params:
        outputfolder = base_path + "/4_R_HD"
    conda:
        "envs/R.yaml"
    threads: thread_ultrasmall
    shell:
        '''
        Rscript R/HD_R.R --input {input} --output {params.outputfolder} 2> {log.out} 1> {log.err}
        '''


rule Cas10_GGDD_hmm_maker:
    '''
    1. Extracts the X AA of sequence around GGDDs in Cas10s.
    2. Aligns the Cas10s using Muscle
    3. Creates a new HD-HMM profile from this alignment

    As of 23.8.2023, running the script two times: once for GGDD and once for GGDE.
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        cas10_sequences = rules.Cas10_concatenate.output,
    output:
        hmm = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.hmm",
        msa = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.msa",
        faa = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.faa",
        hmm_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.hmm",
        msa_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.msa",
        faa_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.faa",
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    shell:
        '''
        python scripts/Cas10_GGDD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm} --msa {output.msa} --faa {output.faa} --motif "GGDD"
        python scripts/Cas10_GGDD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm_GGDE} --msa {output.msa_GGDE} --faa {output.faa_GGDE} --motif "GGDE"
        touch {output.msa}
        touch {output.faa}
        '''

rule Cas10_GGDD_hmmer:
    '''
    Using the GGDD-HMM profiles made in rule Cas10_GGDD_hmm_maker, searches for GGDD domains in all Cas10s
    '''
    input:
        hmm_db_GGDD = rules.Cas10_GGDD_hmm_maker.output.hmm,
        hmm_db_GGDE = rules.Cas10_GGDD_hmm_maker.output.hmm_GGDE,
        cas10_sequences = rules.Cas10_concatenate.output
    output:
        raw_table_GGDD = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits.out",
        raw_table_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM_hits.out"
    params:
        out_GGDD = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_temp.out",
        out_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM_hits_temp.out",
        rows_GGDD = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_rows.out",
        rows_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM_hits_rows.out",
        rows1 = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_rows1.out",
        E = "1e-02"
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    shell:
        '''
        hmmscan --tblout {params.out_GGDD} --cpu 8 -E {params.E} {input.hmm_db_GGDD} {input.cas10_sequences}
        grep -v "#" {params.out_GGDD} > {params.rows_GGDD}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table_GGDD}
        cat {params.rows_GGDD} >> {output.raw_table_GGDD}

        hmmscan --tblout {params.out_GGDE} --cpu 8 -E {params.E} {input.hmm_db_GGDE} {input.cas10_sequences}
        grep -v "#" {params.out_GGDE} > {params.rows_GGDE}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table_GGDE}
        cat {params.rows_GGDE} >> {output.raw_table_GGDE}
        '''

rule merge_GGDD_hmm_with_info:
    '''
    Merges the GGDD-HMM hits with the type III info table
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        GGDD_hmm_hits = rules.Cas10_GGDD_hmmer.output.raw_table_GGDD,
        GGDE_hmm_hits = rules.Cas10_GGDD_hmmer.output.raw_table_GGDE,
        cas10s_faa = rules.Cas10_concatenate.output
    output:
        merged_table = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_merged.tsv"
    params:
        output_folder = base_path + "/Cas10_GGDD_hmm_profiles"
    threads: thread_hogger
    run:
        import pandas as pd
        from Bio import SeqIO

        info = pd.read_csv(str(input.info_table), sep = "\t")
        hmm_hits_GGDD = pd.read_csv(str(input.GGDD_hmm_hits), sep = '\s+')
        hmm_hits_GGDE = pd.read_csv(str(input.GGDE_hmm_hits), sep = '\s+')

        #reduce the hmm_hits to only query_name and E-value_best_domain
        hmm_hits_GGDD = hmm_hits_GGDD[["query_name", "E-value_best_domain"]]
        hmm_hits_GGDE = hmm_hits_GGDE[["query_name", "E-value_best_domain"]]

        #rename the columns
        hmm_hits_GGDD.columns = ["Locus", "GGDD_E-value"]
        hmm_hits_GGDE.columns = ["Locus", "GGDE_E-value"]
        #add column "GGDD_hmm_boolean" and set it to True
        hmm_hits_GGDD["GGDD_hmm_boolean"] = True
        hmm_hits_GGDE["GGDE_hmm_boolean"] = True

        #merge the info table with the hmm_hits_GGDD table
        merged = pd.merge(info, hmm_hits_GGDD, on = "Locus", how = "left")
        merged = pd.merge(merged, hmm_hits_GGDE, on = "Locus", how = "left")

        #if some rows were left without a hit, set the GGDD_hmm_boolean to False
        merged["GGDD_hmm_boolean"] = merged["GGDD_hmm_boolean"].fillna(False)
        merged["GGDE_hmm_boolean"] = merged["GGDE_hmm_boolean"].fillna(False)


        #read in the Cas10 sequences using biopython
        cas10s = SeqIO.to_dict(SeqIO.parse(str(input.cas10s_faa), "fasta"))

        all_cyclase_literals = ["SGDD", "AGDD", "GGED", "GGDD", "GGED", "GGDE", "EGDD", "GEDD", "DGDD", "AGDE", "KGDD", "AGDE"] #updated list of possible cyclase motifs

        #create a dictionary of the literals and add count as value
        all_cyclase_literals_dic = {literal: 0 for literal in all_cyclase_literals}

        #create a column cyclase_literal_sequence in the merged table
        merged["cyclase_literal_sequence"] = ""

        #the merged contains a column "cyclase_literal". This column is true if the sequence contains one of the motifs in all_cyclase_literals. Go through each cas10 sequence
        for index, row in merged.iterrows():
            locus = row["Locus"]
            try:
                cas10 = cas10s[locus].seq
                #check if the sequence contains any of the motifs
                for motif in all_cyclase_literals:
                    if motif in cas10:
                        merged.at[index, "cyclase_literal"] = True
                        print("Cyclase found in " + locus + " with motif " + motif)
                        merged.at[index, "cyclase_literal_sequence"] = motif
                        break
                else:
                    print("No cyclase found for " + locus)
                    merged.at[index, "cyclase_literal"] = False
            except:
                print("No Cas10 found in fasta for " + locus)

        #create new column cyclase which is true only if column GGDD_hmm_boolean or column GGDE_hmm_boolean is true and cyclase_literal is true
        merged["cyclase"] = merged["cyclase_literal"] & (merged["GGDD_hmm_boolean"] | merged["GGDE_hmm_boolean"])

        #in the all_cyclase_literals_dic, add the count of each motif
        for index, row in merged.iterrows():
            motif = row["cyclase_literal_sequence"]
            if motif in all_cyclase_literals_dic:
                all_cyclase_literals_dic[motif] += 1

        #create a new dataframe from the dictionary
        cyclase_literal_df = pd.DataFrame(list(all_cyclase_literals_dic.items()), columns = ["motif", "count"])

        #count percentages of each motif
        cyclase_literal_df["percentage"] = round(cyclase_literal_df["count"] / len(merged) * 100, 1)

        #save the dataframe to a file
        cyclase_literal_df.to_csv(params.output_folder + "/cyclase_literal_count.tsv", sep = "\t", index = False)

        #from the final file, remove all columns but Locus, GGDD_E-value and GGDD_hmm_boolean and cyclase_literal_sequence
        merged = merged[["Locus", "GGDD_E-value", "GGDD_hmm_boolean", "GGDE_E-value", "GGDE_hmm_boolean", "cyclase", "cyclase_literal", "cyclase_literal_sequence"]]

        #save the merged table
        merged.to_csv(str(output.merged_table), sep = "\t", index = False)


rule combine_GGDD_HMM_to_mastertable:
    '''
    Merges the HD and GGDD-HMM data with the master table. Also merges domains info.
    Note that this is not the final mastertable (see rule mastercombiner)
    '''
    input:
        info = rules.concatenate_type_iii_info.output,
        HD_hmm_hits = rules.merge_HD_hmm_with_info.output.merged_table,
        GGDD_hmm_hits = rules.merge_GGDD_hmm_with_info.output.merged_table,
        domains = rules.annotate_bacteria_and_archaea_domains.output,
    output:
        final_info_table = base_path + "/mastertable.tsv"
    threads: thread_ultrasmall
    run:
        import pandas as pd
        info = pd.read_csv(str(input.info), sep = "\t")
        domains = pd.read_csv(str(input.domains), sep = "\t")
        hmm_hits_HD = pd.read_csv(str(input.HD_hmm_hits), sep = "\t")
        hmm_hits_GGDD = pd.read_csv(str(input.GGDD_hmm_hits), sep = "\t")
        print(hmm_hits_GGDD)
        merged = pd.merge(info, hmm_hits_HD, on = "Locus", how = "left")
        merged = pd.merge(merged, hmm_hits_GGDD, on = "Locus", how = "left")
        merged = pd.merge(merged, domains, left_on = "Sample", right_on = "sample", how = "left")
        merged.to_csv(str(output.final_info_table), sep = "\t", index = False)


#a rule that uses the output file cora_type_iii_info from the rule above filter samples with CorA and further divide them into CRISPR-Cas subtypes
rule cora_plot_extractor:
    '''
    Gets cctyper generated plots for each sample with a CorA and a type III CRISPR-Cas system.
    '''
    input: 
        info = rules.concatenate_type_iii_info.output[0]
    output:
        done = base_path + "/xtra1_cora_iii_loci_plots/done.done",
    threads: thread_small
    run:
        import pandas as pd
        import shutil
        #using pandas, get the info from the cora_type_iii_info file and filter out samples that do not have a CorA
        cora_type_iii_info = pd.read_csv(input.info, sep = "\t")
        cora_type_iii_info = cora_type_iii_info[cora_type_iii_info["CorA"] == True]

        #create output folder if it does not exist
        if not os.path.exists(base_path + "/xtra1_cora_iii_loci_plots"):
            os.makedirs(base_path + "/xtra1_cora_iii_loci_plots")

        #For each sample in the pandas dataframe, extract the corresponding plot from the cctyper folder (path is base_path + "/07_cctyper/{j}/plot.png) and copy it to the output folder in a CRISPR-Cas subtype subfolder (e.g. III-A or III-B) depending on the Subtype column in the cora_type_iii_info file
        for index, row in cora_type_iii_info.iterrows():
            sample = row["Sample"]
            subtype = row["Subtype"]
            print(sample + "," + subtype)
            #check if subtype does not contain the substring "Hybrid" (this is because some samples have a hybrid subtype, e.g. III-A/B, and we don't want to create a III-A/B folder)
            if "Hybrid" not in subtype:
                #if the subtype folder does not exist, create it
                if not os.path.exists(base_path + "/xtra1_cora_iii_loci_plots/" + subtype):
                    os.makedirs(base_path + "/xtra1_cora_iii_loci_plots/" + subtype)
                shutil.copyfile(base_path + "/07_cctyper/" + sample + "/plot.png", base_path + "/xtra1_cora_iii_loci_plots/" + subtype + "/" + sample + "_plot.png")

        #create the done.done file to indicate that the rule has finished running
        open(output.done, "w").close()


rule CorA_concatenate:
    '''
    This version works on the cATyper outputs instead of CCTyper outputs
    '''
    input:
        coras = aggregate_CorA_sequences_catyper,
        type_iii_wildcarder_finished = rules.concatenate_type_iii_info.output
    output: base_path + "/09_crispr_iii_CorA/CorAs.faa"
    threads: thread_ultrasmall
    shell:
        '''
        cat {input.coras} > {output}
        '''

rule CorA_cluster:
    '''

    '''
    input: rules.CorA_concatenate.output
    output:
        proteins = base_path + "/15_CorA_cluster/CorA_cluster.faa",
        clusterinfo = base_path + "/15_CorA_cluster/CorA_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/15_CorA_cluster/logs/CorA_align.out",
        err = base_path + "/15_CorA_cluster/logs/CorA_align.err"
    threads: thread_hogger
    shell:
        '''
        cd-hit -i {input} -o {output.proteins} -c 0.90 -n 5 -d 0 -M 16000 -T {threads}
        '''


rule CorA_align:
    input: rules.CorA_cluster.output.proteins
    output: base_path + "/16_CorA_align/CorA_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/16_CorA_align/logs/CorA_align.out",
        err = base_path + "/16_CorA_align/logs/CorA_align.err"
    threads: thread_hogger
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule CorA_align_unclustered:
    input: rules.CorA_concatenate.output
    output: base_path + "/16_CorA_align_unclustered/CorA_alignment_unclustered.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/16_CorA_align/logs/CorA_align.out",
        err = base_path + "/16_CorA_align/logs/CorA_align.err"
    threads: thread_hogger
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule CorA_tree:
    '''

    '''
    input: rules.CorA_align.output
    output: base_path + "/17_CorA_tree/CorA_tree.txt",
    conda: "envs/trees.yaml"
    threads: thread_hogger
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule CorA_tree_unclustered:
    '''

    '''
    input: rules.CorA_align_unclustered.output
    output: base_path + "/17_CorA_tree_unclustered/CorA_tree_unclustered.txt",
    conda: "envs/trees.yaml"
    threads: thread_hogger
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule Cas10_cluster:
    '''
    '''
    input: 
        cas10 = rules.Cas10_concatenate.output
    output:
        proteins = base_path + "/10_Cas10_cluster/cas10_cluster.faa",
        clusterinfo = base_path + "/10_Cas10_cluster/cas10_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/10_Cas10_cluster/logs/cas10_align.out",
        err = base_path + "/10_Cas10_cluster/logs/cas10_align.err"
    threads: thread_hogger
    shell:
        '''
        echo {protein_clustering}
        if [ {protein_clustering} = "True" ]; then
            cd-hit -i {input.cas10} -o {output.proteins} -c 0.99 -n 5 -d 0 -M 16000 -T {threads}
        else
            cp {input.cas10} {output.proteins}
            touch {output.clusterinfo}
        fi
        '''


rule Cas10_align:
    '''
    
    '''
    input: rules.Cas10_cluster.output.proteins
    output: base_path + "/11_Cas10_align/cas10_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/11_Cas10_align/logs/cas10_align.out",
        err = base_path + "/11_Cas10_align/logs/cas10_align.err"
    threads: thread_hogger
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule Cas10_tree:
    '''

    '''
    input: rules.Cas10_align.output
    output: base_path + "/12_Cas10_tree/cas10_tree.txt",
    conda: "envs/trees.yaml"
    threads: thread_hogger
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule validated_new_effectors:
    '''
    This rule is updated as potential new effectors from group4 emerge.
    Creates tables for each locus showing whether a validated effector is
    found in that locus. Does not interfere with running the group4 analysis
    and the validated effectors found here will also be listed in group4 proteins,
    as that is where they were originally found.
    '''
    input:
        proteins = rules.cATyper_hmm_search.output.contig_proteins
    output:
        #validated_new_effectors = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_orig.tsv",
        temp_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_temp.tsv",
        hmm_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        validated_effectors_hmm_db = validated_effectors_hmm_db,
        temp_hmm = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/60_validated_new_effectors/logs/{c}/{c}_validated_new_effectors.out",
        err = base_path + "/60_validated_new_effectors/logs/{c}/{c}_validated_new_effectors.err"
    threads: thread_ultrasmall
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.validated_effectors_hmm_db} {input.proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

rule concatenate_validate_new_effectors_hmm:
    input: aggregate_validated_new_effectors_hmm
    output: base_path + "/60_validated_new_effectors/validated_new_effectors_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule validated_new_effectors_analysis:
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv",
    output:
        catyper = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/61_validated_new_effectors_analysis/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = validated_effectors_folder + "/profiles",
        tmhmm_model_path = TM_path + "/TMHMM2.0.model",
    conda: "envs/groupChar.yaml"
    log:
        out = base_path + "/61_validated_new_effectors_analysis/logs/{c}.out",
        err = base_path + "/61_validated_new_effectors_analysis/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        echo "Running validated effector analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "validated_new_effector" --effector_plot_data {output.plottable_effector_positions} --tmhmm_model_path {params.tmhmm_model_path} 2> {log.err} 1> {log.out}
        '''

rule concatenate_validated_new_effectors_analysis:
    input: aggregate_validated_new_effectors_analysis
    output: base_path + "/61_validated_new_effectors_analysis/validated_new_effectors_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """


rule concatenate_validated_new_effectors_scores:
    '''
    Concatenates info on the effector hits from cATyper analysis
    '''
    input:
        protein_to_effector = aggregate_validated_new_effectors_pte,
        effector_to_protein = aggregate_validated_new_effectors_etp
    output:
        protein_to_effector_concatenated = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_to_protein.tsv",
    threads: thread_ultrasmall
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule analyse_validated_new_effectors_scores:
    '''
    The algorithm that analyses how well each hmm profile performs in validated_new_effectors.
    '''
    input:
        pte = rules.concatenate_validated_new_effectors_scores.output.protein_to_effector_concatenated,
        etp = rules.concatenate_validated_new_effectors_scores.output.effector_to_protein_concatenated,
    output:
        effector_scores = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_scores.tsv",
        effector_scores_summary = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_scores_summary.tsv",
        #effector_scores_plot1 = base_path + "/11_cATyper_analysis/cATyper_effector_scores_plot1.png",
    params:
        outdir = base_path + "/61_validated_new_effectors_analysis"
    threads: thread_ultrasmall
    conda: "envs/analysis.yaml"
    shell:
        '''
        python scripts/catyper_effector_scores.py --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        touch {output.effector_scores}
        touch {output.effector_scores_summary}
        '''

rule heatmap_known_validated_effectors:
    input:
        etp_validated = rules.concatenate_validated_new_effectors_scores.output.effector_to_protein_concatenated,
        etp_known = rules.concatenate_cATyper_analysis_effector_scores.output.effector_to_protein_concatenated
    output:
        effectors_combined = base_path + "/70_validated_and_known_effectors_heatmap/etp_combined.tsv",
        effector_scores_summary = base_path + "/70_validated_and_known_effectors_heatmap/effectors.tsv",
    params:
        outdir = base_path + "/70_validated_and_known_effectors_heatmap",
        dummyout = base_path + "/70_validated_and_known_effectors_heatmap/dummyout.tsv"
    threads: thread_ultrasmall
    conda: "envs/analysis.yaml"
    shell:
        '''
        cat {input.etp_validated} {input.etp_known} > {output.effectors_combined}
        python scripts/catyper_effector_scores.py --etp {output.effectors_combined} --output {params.dummyout} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        '''


#### RING NUCLEASE RULES ####

rule ring_nucleases:
    '''
    Reuses code from validated new effectors to find ring nucleases.
    Runs hmmscan against all loci using the custom ring nuclease database as, well, the database.
    '''
    input:
        proteins = rules.cATyper_hmm_search.output.contig_proteins
    output:
        temp_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_temp.tsv",
        hmm_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        rn_db = ring_nuclease_db,
        temp_hmm = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases.temp",
        evalue = "1e-8"
    log:
        out = base_path + "/70_ring_nucleases/logs/{c}/{c}_ring_nucleases.out",
        err = base_path + "/70_ring_nucleases/logs/{c}/{c}_ring_nucleases.err"
    threads: thread_ultrasmall
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.rn_db} {input.proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

rule concatenate_validate_ring_nucleases_hmm:
    input: aggregate_ring_nucleases_hmm
    output: base_path + "/70_ring_nucleases/ring_nucleases_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule ring_nucleases_analysis:
    '''
    Generates locus-specific information on the presence of ring nucleases.
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv",
    output:
        catyper = base_path + "/71_ring_nucleases_analysis/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/71_ring_nucleases_analysis/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/71_ring_nucleases_analysis/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/71_ring_nucleases_analysis/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/71_ring_nucleases_analysis/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/71_ring_nucleases_analysis/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = ring_nuclease_folder + "/profiles",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/71_ring_nucleases_analysis/logs/{c}.out",
        err = base_path + "/71_ring_nucleases_analysis/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        echo "Running validated effector analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "ring_nucleases" --effector_plot_data {output.plottable_effector_positions} --ring_nuclease True 2> {log.err} 1> {log.out}
        '''

rule ring_nucleases_analysis_no_length_limit:
    '''
    Generates locus-specific information on the presence of ring nucleases.
    This version of the rule has no length limitations, allowing the search for fusions of RNs with other proteins.
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv",
    output:
        catyper = base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/71_ring_nucleases_analysis_no_length_lim/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = ring_nuclease_folder + "/profiles",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/71_ring_nucleases_analysis_no_length_lim/logs/{c}.out",
        err = base_path + "/71_ring_nucleases_analysis_no_length_lim/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        echo "Running validated effector analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "ring_nucleases" --effector_plot_data {output.plottable_effector_positions} --ring_nuclease True --no_length_lim True 2> {log.err} 1> {log.out}
        '''

rule concatenate_ring_nucleases_analysis:
    input: aggregate_ring_nucleases_analysis
    output: base_path + "/71_ring_nucleases_analysis/ring_nucleases_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule concatenate_ring_nucleases_analysis_no_length_lim:
    input: aggregate_ring_nucleases_analysis_no_length_lim
    output: base_path + "/71_ring_nucleases_analysis_no_length_lim/ring_nucleases_all_hmm_no_length_lim.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule concatenate_ring_nucleases_scores:
    '''
    Concatenates info on the ring nuclease hits from the ring nuclease analysis
    '''
    input:
        protein_to_effector = aggregate_ring_nucleases_pte,
        effector_to_protein = aggregate_ring_nucleases_etp,
        ring_nuclease_concat = rules.concatenate_ring_nucleases_analysis.output
    output:
        protein_to_effector_concatenated = base_path + "/71_ring_nucleases_analysis/ring_nucleases_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/71_ring_nucleases_analysis/ring_nucleases_effectors_to_protein.tsv",
    threads: thread_small
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule concatenate_ring_nucleases_scores_no_length_lim:
    '''
    Concatenates info on the ring nuclease hits from the ring nuclease analysis.
    No length limitation version.
    '''
    input:
        protein_to_effector = aggregate_ring_nucleases_pte_no_length_lim,
        effector_to_protein = aggregate_ring_nucleases_etp_no_length_lim,
        ring_nuclease_concat = rules.concatenate_ring_nucleases_analysis_no_length_lim.output
    output:
        protein_to_effector_concatenated = base_path + "/71_ring_nucleases_analysis_no_length_lim/ring_nucleases_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/71_ring_nucleases_analysis_no_length_lim/ring_nucleases_effectors_to_protein.tsv",
    threads: thread_small
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule analyse_ring_nucleases_scores:
    '''
    The algorithm that analyses how well each hmm profile performs in ring_nucleases.
    '''
    input:
        pte = rules.concatenate_ring_nucleases_scores.output.protein_to_effector_concatenated,
        etp = rules.concatenate_ring_nucleases_scores.output.effector_to_protein_concatenated,
    output:
        effector_scores = base_path + "/71_ring_nucleases_analysis/val_new_effectors_scores.tsv",
        effector_scores_summary = base_path + "/71_ring_nucleases_analysis/val_new_effectors_scores_summary.tsv",
        #effector_scores_plot1 = base_path + "/11_cATyper_analysis/cATyper_effector_scores_plot1.png",
    params:
        outdir = base_path + "/71_ring_nucleases_analysis"
    threads: thread_small
    conda: "envs/analysis.yaml"
    shell:
        '''
        python scripts/catyper_effector_scores.py --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        touch {output.effector_scores}
        touch {output.effector_scores_summary}
        '''


rule heatmap_ring_nucleases:
    input:
        ring_nuclease_etp = rules.concatenate_ring_nucleases_scores.output.effector_to_protein_concatenated
    output:
        rns_combined = base_path + "/72_ring_nucleases_heatmap/etrn_combined.tsv", # etrn = effector to ring nuclease
        rn_scores_summary = base_path + "/72_ring_nucleases_heatmap/rns.tsv",
    params:
        outdir = base_path + "/72_ring_nucleases_heatmap",
        dummyout = base_path + "/72_ring_nucleases_heatmap/dummyout.tsv"
    threads: thread_small
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/72_ring_nucleases_heatmap/logs/heatmap.out",
        err = base_path + "/72_ring_nucleases_heatmap/logs/heatmap.err"
    shell:
        '''
        python scripts/catyper_effector_scores.py --etp {input.ring_nuclease_etp} --output {params.dummyout} --output_summary {output.rn_scores_summary} --outdir {params.outdir} 2> {log.err} 1> {log.out}
        '''

rule ring_nucleases_fusions:
    '''
    Looks at the original ring nuclease hits and finds out if a protein has multiple credible ring nuclease annotations
    to indicate a fusion between two ring nucleases.
    '''
    input:
        hmm_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv",
    output:
        ring_nuclease_fusions = base_path + "/73_ring_nuclease_fusions/{c}/{c}_ring_nuclease_fusions.tsv",
    params:
        outdir = base_path + "/73_ring_nuclease_fusions/{c}",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/73_ring_nuclease_fusions/logs/{c}.out",
        err = base_path + "/73_ring_nuclease_fusions/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/ringnucleasefusionfinder.py --locus {wildcards.c} --output_folder {params.outdir} --hmm_rows {input.hmm_rows} 2> {log.err} 1> {log.out}
        '''

rule concatenate_ring_nuclease_fusions:
    input: aggregate_ring_nuclease_fusions
    output: base_path + "/73_ring_nuclease_fusions/ring_nuclease_fusions_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """


rule ring_fusions:
    '''
    Takes protein_to_effector data from all three search modules (known effectors, new validated effectors and ring nucleases)
    and by comparing the proteins listed in ring nucleases, creates a table
    showing to which effectors the ring nucleases are fused to
    '''
    input:
        known_effectors = rules.concatenate_cATyper_analysis_effector_scores.output.protein_to_effector_concatenated,
        validated_new_effectors = rules.concatenate_validated_new_effectors_scores.output.protein_to_effector_concatenated,
        ring_nucleases = rules.concatenate_ring_nucleases_scores_no_length_lim.output.protein_to_effector_concatenated,
    output:
        ring_fusions = base_path + "/100_ring_fusions/ring_fusions.tsv"
    threads: thread_ultrasmall
    log:
        out = base_path + "/100_ring_fusions/logs/ring_fusions.out",
        err = base_path + "/100_ring_fusions/logs/ring_fusions.err"
    shell:
        '''
        python scripts/ring_fusions.py --known_effectors {input.known_effectors} --validated_new_effectors {input.validated_new_effectors} --ring_nucleases {input.ring_nucleases} --output {output.ring_fusions} 2> {log.err} 1> {log.out}
        '''

rule ring_nuclease_cas10_fusions:
    '''
    Runs Cas10 sequences as query against the ring nucleases database to detect any fusions
    '''
    input:
        cas10 = rules.Cas10_cluster.output.proteins
    output:
        temp_rows = base_path + "/101_ring_fusions_cas10/ring_nuclease_cas10_temp.tsv",
        hmm_rows = base_path + "/101_ring_fusions_cas10/ring_nuclease_cas10_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        ring_nuclease_hmm_db = ring_nuclease_db,
        temp_hmm = base_path + "/101_ring_fusions_cas10/ring_nucleases_cas10.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/101_ring_fusions_cas10/logs/ring_fusions_cas10.out",
        err = base_path + "/101_ring_fusions_cas10/logs/ring_fusions_cas10.err"
    threads: thread_hogger
    shell:
        '''        
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.ring_nuclease_hmm_db} {input.cas10}
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        cat {output.temp_rows} >> {output.hmm_rows}
        '''

rule validated_effectors_cas10_fusions:
    '''
    Runs Cas10 sequences as query against the validated effectors database to detect any fusions
    '''
    input:
        cas10 = rules.Cas10_cluster.output.proteins
    output:
        temp_rows = base_path + "/101_validated_effectors_cas10/validated_effectors_cas10_temp.tsv",
        hmm_rows = base_path + "/101_validated_effectors_cas10/validated_effectors_cas10_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        validated_effectors_hmm_db = validated_effectors_hmm_db,
        temp_hmm = base_path + "/101_validated_effectors_cas10/ring_nucleases_cas10.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/101_validated_effectors_cas10/logs/validated_effectors_cas10.out",
        err = base_path + "/101_validated_effectors_cas10/logs/validated_effectors_cas10.err"
    threads: thread_hogger
    shell:
        '''        
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.validated_effectors_hmm_db} {input.cas10}
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        cat {output.temp_rows} >> {output.hmm_rows}
        '''

rule known_effectors_cas10_fusions:
    '''
    Runs Cas10 sequences as query against the known effectors database to detect any fusions
    '''
    input:
        cas10 = rules.Cas10_cluster.output.proteins
    output:
        temp_rows = base_path + "/101_known_effectors_cas10/known_effectors_cas10_temp.tsv",
        hmm_rows = base_path + "/101_known_effectors_cas10/known_effectors_cas10_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        known_effectors_db = hmm_database_folder + "/" + hmm_database_file,
        temp_hmm = base_path + "/101_known_effectors_cas10/ring_nucleases_cas10.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/101_known_effectors_cas10/logs/known_effectors_cas10.out",
        err = base_path + "/101_known_effectors_cas10/logs/known_effectors_cas10.err"
    threads: thread_hogger
    shell:
        '''        
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.known_effectors_db} {input.cas10}
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        cat {output.temp_rows} >> {output.hmm_rows}
        '''

rule mastercombiner:
    '''
    Added 4.9.2023 to incorporate merging and data handling previously done in R in here
    '''
    input:
        concat_taxInfo = rules.concat_taxInfo.output,
        catyper = rules.concatenate_cATyper_analysis.output,
        type_iii_info = rules.combine_GGDD_HMM_to_mastertable.output.final_info_table,
        temperature_data = temperature_data,
        validated_new_effectors = rules.concatenate_validated_new_effectors_analysis.output,
        ring_fusions = rules.ring_fusions.output.ring_fusions,
        ring_nucleases = rules.concatenate_ring_nucleases_analysis.output
    output:
        final_info_table = base_path + "/mastertable_v2.tsv",
        group4_IDs = base_path + "/group4/group4_IDs.txt"
    threads: thread_small
    run:
        import pandas as pd
        mastertable = pd.read_csv(str(input.type_iii_info), sep = "\t")
        taxInfo = pd.read_csv(str(input.concat_taxInfo), sep = "\t")
        catyper = pd.read_csv(str(input.catyper), sep = "\t")
        validated_new_effectors_df = pd.read_csv(str(input.validated_new_effectors), sep = "\t")
        ring_nucleases_df = pd.read_csv(str(input.ring_nucleases), sep = "\t")

        #### RING FUSION PREP ####
        ring_fusions = pd.read_csv(str(input.ring_fusions), sep = "\t")
        #in ring fusions only retain columns locus, ring_nuclease, fusion_components, fusion_protein
        ring_fusions = ring_fusions[["locus", "ring_nuclease", "fusion_components", "fusion_protein"]]

        #### KNOWN EFFECTORS PREP ####
        #in catyper, if any of the following columns are True, set has_known_effector column to True: ca3	ca4	ca5	ca6	sam-amp
        catyper["has_known_effector"] = False
        catyper.loc[catyper[["ca3", "ca4", "ca5", "ca6", "sam-amp"]].any(axis = 1), "has_known_effector"] = True

        #### VALIDATED NEW EFFECTOR PREP ####
#       validated_new_effectors_df = validated_new_effectors_df[["can3", "tirsaved", "cam2", "mem-01", "mem2", "locus"]] #drop cOA columns
        validated_new_effectors_df["has_validated_new_effector"] = False
        #if any column except locus is True, set has_validated_new_effector to True
        validated_new_effectors_df.loc[validated_new_effectors_df[["can3", "tirsaved", "cam3", "cam2"]].any(axis = 1), "has_validated_new_effector"] = True

        #rename column no_effectors to no_validated_new_effectors
        validated_new_effectors_df = validated_new_effectors_df.rename(columns = {"no_effectors": "no_validated_new_effectors"})

        #rename columns ca3	ca4	ca5	ca6 to NE_ca3 NE_ca4 NE_ca5 NE_ca6
        validated_new_effectors_df = validated_new_effectors_df.rename(columns = {"ca3": "NE_ca3", "ca4": "NE_ca4", "ca5": "NE_ca5", "ca6": "NE_ca6", "sam-amp": "NE_sam-amp"})

        #### RING NUCLEASE PREP ####
        ring_nucleases_df["has_ring_nuclease"] = False
        #if any of the following are true, set has_ring_nuclease to True
        ring_nucleases_df.loc[ring_nucleases_df[list_of_RNs].any(axis = 1), "has_ring_nuclease"] = True
        
        #drop ca4, ca4, ca6 sam-amp and unk from ring_nucleases
        ring_nucleases_df = ring_nucleases_df.drop(columns = ["ca3", "ca4", "ca5", "ca6", "sam-amp", "unk", "mem", "no_effectors"])

        #### MERGE ####
        mastertable = pd.merge(mastertable, taxInfo, on = "Sample", how = "left")
        mastertable = pd.merge(mastertable, catyper, left_on = "Locus", right_on = "locus", how = "left")
        mastertable = pd.merge(mastertable, validated_new_effectors_df, left_on = "Locus", right_on = "locus", how = "left", suffixes = ("", "new_eff"))
        mastertable = pd.merge(mastertable, ring_nucleases_df, left_on = "Locus", right_on = "locus", how = "left", suffixes = ("", "rn"))
        mastertable = pd.merge(mastertable, ring_fusions, left_on = "Locus", right_on = "locus", how = "left", suffixes = ("", "rn_fusions"))

        #create columns con_ca3, con_ca4, con_ca5, con_ca6, con_sam_amp. Con stands for consensus
        mastertable["con_ca3"] = False
        mastertable["con_ca4"] = False
        mastertable["con_ca5"] = False
        mastertable["con_ca6"] = False
        mastertable["con_sam-amp"] = False
        #print columsn from mastertable
        #print(mastertable.columns)
        #if ca3 or NE_ca3 is True, set con_ca3 to True. Same for all signal molecules
        mastertable.loc[(mastertable["ca3"] == True) | (mastertable["NE_ca3"] == True), "con_ca3"] = True
        mastertable.loc[(mastertable["ca4"] == True) | (mastertable["NE_ca4"] == True), "con_ca4"] = True
        mastertable.loc[(mastertable["ca5"] == True) | (mastertable["NE_ca5"] == True), "con_ca5"] = True
        mastertable.loc[(mastertable["ca6"] == True) | (mastertable["NE_ca6"] == True), "con_ca6"] = True
        mastertable.loc[(mastertable["sam-amp"] == True) |(mastertable["NE_sam-amp"]), "con_sam-amp"] = True

        mastertable["multiple_signals"] = False
        #if more than one of the con columns is True, set multiple_signals to True
        mastertable.loc[(mastertable[["con_ca3", "con_ca4", "con_ca5", "con_ca6", "con_sam-amp"]].sum(axis = 1) > 1), "multiple_signals"] = True

        mastertable["effector_count_known_new_sum"] = 0
        #sum the number of effectors in the known and validated new effectors columns. Effectors to sum are can3, tirsaved, cam2, cam3, cam2, saved-chat, nucc, calpl, cami1, can1-2, csx1, csx23 csm6-ca6, cora
        mastertable["effector_count_known_new_sum"] = mastertable[["can3", "tirsaved", "cam3", "cam2", "saved-chat", "nucc", "calpl", "cami1", "can1-2", "csx1", "csx23", "csm6-ca6", "cora", "cam1"]].sum(axis = 1)

        #for mastertable rows in which fusion_protein is NaN, make it False
        mastertable["fusion_protein"] = mastertable["fusion_protein"].fillna(False)

        #open the temperature data file
        temperature_data = pd.read_csv(str(input.temperature_data), sep = ",")
        #group by column "genus" and calculate mean temperature for each genus
        mean_temperatures = temperature_data.groupby("genus")["Topt_ave"].mean().dropna()
        print(mean_temperatures)
        #merge mean_temperature with mastertable
        mastertable = pd.merge(mastertable, mean_temperatures, on = "genus", how = "left")
        #rename topt_ave to "temperature"
        mastertable = mastertable.rename(columns = {"Topt_ave": "mean_temperature"})

        #Tidy the mastertable
        #drop the following columns (list)
        droppable_columns = ["sample_x", "Unnamed: 0", "val_x", "locus_x", "sample_y", "val_y", "locus_y", "Unnamed: 0_y", "mem_y", "Unnamed: 0_x"]

        #drop the columns in a for loop with try/catch for each column in case the mastertable changes in the future
        for column in droppable_columns:
            try:
                mastertable = mastertable.drop(columns = [column])
            except:
                pass

        #create a new column multilocus_sample. This a boolean value that reports whether the same sample has other loci
        mastertable["multilocus_sample"] = False

        #create a new column multisignal_sample. This is a boolean value that reports whether the same sample has multiple signals (ca3, ca4, ca6 or sam-amp)
        mastertable["multisignal_sample"] = False

        #create a new column all_signals_sample. This is a string that contains comma separated list of all signals present in the sample
        mastertable["all_signals_sample"] = ""

        #for each sample, go through loci and list which signal molecules are true and update the multisingal_sample and all_signals_sample columns accordingly
        for sample in mastertable["Sample"].unique():
            #get all loci for the sample
            sample_loci = mastertable[mastertable["Sample"] == sample]["Locus"].unique()
            #if the number of loci is greater than 1, set multilocus_sample to True
            if len(sample_loci) > 1:
                mastertable.loc[mastertable["Sample"] == sample, "multilocus_sample"] = True

            #get all signals for the sample
            sample_signals = mastertable[mastertable["Sample"] == sample][["con_ca3", "con_ca4", "con_ca5", "con_ca6", "con_sam-amp"]].sum(axis = 1)
            #if the number of signals is greater than 1, set multisignal_sample to True
            if sample_signals.sum() > 1:
                mastertable.loc[mastertable["Sample"] == sample, "multisignal_sample"] = True

            if sample_signals.sum() > 0:
                mastertable.loc[mastertable["Sample"] == sample, "signal_sample"] = True

            #using columns ca3, ca4, ca6 and sam-amp, create a string of all signals present in the sample by checking if their boolean value is True
            all_signals = ""
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca3"].any() == True:
                all_signals += "ca3,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca4"].any() == True:
                all_signals += "ca4,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca6"].any() == True:
                all_signals += "ca6,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_sam-amp"].any() == True:
                all_signals += "sam-amp,"
            #remove the last comma from the string
            all_signals = all_signals[:-1]
            #update the all_signals_sample column with the string
            mastertable.loc[mastertable["Sample"] == sample, "all_signals_sample"] = all_signals

        #create column effector_elsewhere_in_sample
        mastertable["effector_elsewhere_in_sample"] = False
        mastertable["effector_elsewhere_in_sample_but_not_here"] = False

        #for each locus, check if the sample has another locus that generates a signal and if so, set cas10_effector_elsewhere_in_sample to True
        for locus in mastertable["Locus"].unique():
            #get the sample for the locus
            sample = mastertable[mastertable["Locus"] == locus]["Sample"].unique()[0]
            #get all loci for the sample
            sample_loci = mastertable[mastertable["Sample"] == sample]["Locus"].unique()
            #if any of the loci in the sample have a signal, set cas10_effector_elsewhere_in_sample to True
            if mastertable[mastertable["Locus"].isin(sample_loci)][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == True:
                mastertable.loc[mastertable["Locus"] == locus, "effector_elsewhere_in_sample"] = True

            #if any of the loci in the sample have a signal but this locus does not have a signal, set cas10_effector_elsewhere_in_sample_but_not_here to True
            #first check if current locus has signal
            if mastertable.loc[mastertable["Locus"] == locus][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == False:
                #if it does not have a signal, check if any other loci in the sample have a signal
                if mastertable[mastertable["Locus"].isin(sample_loci)][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == True:
                    mastertable.loc[mastertable["Locus"] == locus, "effector_elsewhere_in_sample_but_not_here"] = True
            

        #subset the mastertable to include samples where (GGDD_hmm_boolean is True or GGDE_hmm_boolean is True) AND HD_hmm_boolean is False
        mastertable["Cas10_HD_coord"] = pd.to_numeric(mastertable["Cas10_HD_coord"], errors='coerce')
        mastertable_group4 = mastertable[
            ((mastertable["GGDD_hmm_boolean"] == True) | (mastertable["GGDE_hmm_boolean"] == True)) #& 
            #(mastertable["HD_hmm_boolean"] == False) & 
            #((mastertable["Cas10_HD_coord"] > 100) | (mastertable["Cas10_HD"] == False)) & 
            #(mastertable["unknown_proteins"] == True)
        ]
        mastertable_group4["Locus"].to_csv(str(output.group4_IDs), sep = "\t", index = False)

        #some loci are wrongly subtyped. Use this dictionary to manually correct them. Key shows locus, value shows corrected subtype
        subtype_corrections = { #NOTE these may change from run to run as crispr locus numbers are arbitrary and stochastic. Any rerun of CCTyper will null these corrections.
            # "GCA_003201765.2_3": "III-D",
            # "GCF_003201765.2_2": "III-D",
            # "GCA_022845955.1_1": "III-D",
            # "GCA_001548075.1_0": "III-D",
             "GCF_014495845.1_1": "III-D",
            # "GCA_000143965.1_0": "III-D",
            # "GCA_000010725.1_0": "III-D",
             "GCF_019900845.1_1": "III-D",
            # "GCF_000227745.2_1": "III-B",
            # "GCA_000970265.1_31": "III-F",
            # "GCF_000022425.1_4": "III-F",
            # "GCF_000022405.1_3": "III-F",
            # "GCF_009602405.1_0": "III-F",
            # "GCA_003201675.2_2": "III-F",
            # "GCF_003201675.2_3": "III-F",
            # "GCF_028472365.1_2": "III-F",
            # "GCA_000338775.1_1": "III-F",
            # "GCA_000011205.1_3": "III-F",
            # "GCA_003967175.1_0": "III-F",
            
            # added 13.8.2024
            "GCF_000022425.1_4": "III-B",
            "GCF_000227745.2_2": "III-B",

            "GCF_003967175.1_2": "III-F",
            "GCF_000022465.1_3": "III-F",
            "GCF_000011205.1_0": "III-F",
            "GCF_000338775.1_1": "III-F",
            "GCF_000340315.1_1": "III-F",
            "GCF_003201675.2_1": "III-F",
            "GCF_012222305.1_3": "III-F",
            "GCF_000022425.1_3": "III-F",
            "GCF_000022405.1_0": "III-F",

            "GCF_003201675.2_3": "III-D",
            "GCF_000143965.1_1": "III-D",
            "GCF_000010725.1_1": "III-D",
            "GCF_963678605.1_1": "III-D",
            "GCF_963677725.1_0": "III-D",

            "GCF_009602405.1_0": "III-F",

            "GCF_902827225.1_0": "III-C",
            "GCF_000017865.1_1": "III-C",
            "GCF_003991135.1_2": "III-C",
            "GCF_000025505.1_0": "III-C",
            "GCF_000970265.1_32": "III-C",

        }

        #correct the subtypes for the given loci (column Locus)
        for locus, subtype in subtype_corrections.items():
            mastertable.loc[mastertable["Locus"] == locus, "Subtype"] = subtype

        #drop any duplicate rows
        mastertable = mastertable.drop_duplicates()
        mastertable.to_csv(str(output.final_info_table), sep = "\t", index = False)


rule node_graph:
    '''
    This rule produces files for Gephi to view effector colocalisation
    using a node graph. It also generates upset plots from the same data.

    For Gephi, it works by first creating a table with all effector combinations (edges) with the headers
    source, target, type and weight. It then looks at each locus in the mastertable
    and upon finding co-occurrence of an effector pair, it adds 1 to the weight column.
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table,
        effector_list = rules.heatmap_known_validated_effectors.output.effectors_combined
    output:
        edges =  base_path + "/80_node_graph/edges_all.tsv",
        nodes = base_path + "/80_node_graph/nodes_all.tsv",
    params:
        outdir = base_path + "/80_node_graph",
        edges_basename = base_path + "/80_node_graph/edges.tsv",
        nodes_basename = base_path + "/80_node_graph/nodes.tsv"
    threads: thread_small
    conda: "envs/effector_nodes.yaml"
    shell:
        '''
        python scripts/effector_nodes_and_other_viz.py --input_mastertable {input.mastertable} --effector_list {input.effector_list} -o {params.outdir} -n {params.nodes_basename} -d {params.edges_basename}
        '''

rule node_graph_RNs:
    '''
    A modification of the rule node_graph. This includes ring nucleases in the node graph and specifically
    draws edges between RNs and effectors.
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table,
        effector_list = rules.heatmap_known_validated_effectors.output.effectors_combined
    output:
        edges =  base_path + "/80_RN_node_graph/edges_all.tsv",
        nodes = base_path + "/80_RN_node_graph/nodes_all.tsv",
    params:
        outdir = base_path + "/80_RN_node_graph",
        edges_basename = base_path + "/80_RN_node_graph/edges.tsv",
        nodes_basename = base_path + "/80_RN_node_graph/nodes.tsv",
        rn_list = list_of_RNs
    threads: thread_small
    conda: "envs/effector_nodes.yaml"
    shell:
        '''
        python scripts/effector_nodes_and_other_viz_RN.py --input_mastertable {input.mastertable} --effector_list {input.effector_list} -o {params.outdir} -n {params.nodes_basename} -d {params.edges_basename}
        '''


rule effector_commonness:
    '''
    This rule counts the number of different effectors across the loci
    '''
    input:
        crispr_loci = rules.mastercombiner.output.final_info_table
    output:
        effector_commonness_tsv = base_path + "/13_effector_commonness/effector_commonness_master.tsv",
        #effector_commonness_plot = base_path + "/13_effector_commonness/effector_commonness.png"
    params:
        base_path = base_path + "/13_effector_commonness/effector_commonness"
    threads: thread_ultrasmall
    shell:
        '''
        python3 scripts/effector_commonness.py --input {input.crispr_loci} --output_basepath {params.base_path}
        '''

rule cctyper_gene_locations:
    '''
    This uses CCTyper outputs to generate a table that shows each cas gene's
    coordinate on the contig and other attributes. All attributes are:
    | protein_id     | start   | end     | effector | locus             | sample           | strand | type
    '''
    input:
        cas_operons_cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        cas_proteins = rules.crispr_locus_proteins.output.crispr_locus_proteins,
    output:
        cctyper_gene_locations_plottable = base_path + "/14_cctyper_gene_locations/{c}/{c}_cctyper_gene_locations.tsv"
    params:
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        this_folder = base_path + "/14_cctyper_gene_locations/{c}",
        outputfolder = base_path + "/14_cctyper_gene_locations/{c}",
        cctyper_folder = base_path + "/07_cctyper"
    log:
        out = base_path + "/14_cctyper_gene_locations/logs/{c}.out",
        err = base_path + "/14_cctyper_gene_locations/logs/{c}.err"
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/cctyper_gene_locations.py --locus_id {wildcards.c} --sample_folder {params.sample_folder} --this_folder {params.this_folder} --outputfolder {params.outputfolder} --cas_operons {input.cas_operons_cctyper} --cctyper_path {params.cctyper_folder} --protein_fasta {input.cas_proteins} 2> {log.err} 1> {log.out}
        '''

rule locus_visualiser:
    '''
    This rule visualises a given range of a gff file.
    Using locus as wildcard.
    GFF comes from genome.
    Coordinates for CRISPR locus (and beyond) comes from 
    Takes CRISPR coordinates from 
    '''
    input:
        #genome_gff = rules.crispr_locus_proteins.output.crispr_locus_gff,
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        cctyper_plottable = rules.cctyper_gene_locations.output.cctyper_gene_locations_plottable,
        plottable_effector_positions = rules.validated_new_effectors_analysis.output.plottable_effector_positions,
        known_effectors = rules.cATyper_analysis.output.plottable_effector_positions,
        ring_nucleases = rules.ring_nucleases_analysis.output.plottable_effector_positions
    output:
        visualisation = base_path + "/90_locus_viz/{c}/{c}_viz.png"
    conda:
        "envs/locus_visualiser.yaml"
    params:
        outdir = base_path + "/90_locus_viz/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        cctyper_folder = base_path + "/07_cctyper",
    log:
        out = base_path + "/90_locus_viz/logs/{c}.out",
        err = base_path + "/90_locus_viz/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/locus_visualiser.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --cctyper_folder {params.cctyper_folder} --validated_effectors {input.plottable_effector_positions} --cctyper_protein_table {input.cctyper_plottable} --known_effector_table {input.known_effectors} --ring_nucleases {input.ring_nucleases}  2> {log.err} 1> {log.out}
        '''

rule concatenate_locus_viz:
    input: aggregate_locus_viz
    output: base_path + "/90_locus_viz/locus_viz.done"
    threads: thread_ultrasmall
    shell:
        '''
        touch {output}
        '''


rule create_html_file:
    '''
    This creates an html file out of the mastertable and adds hyperlinks
    to the loci figures and the associated .tsv files
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table,
    output:
        html_file = base_path + "/type_iii_mastertable.html"
    params:
        viz_dir_full = base_path + "/data",
        viz_dir_relative = "data",
    conda:
        "envs/excel_output.yaml"
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/html_writer.py --mastertable {input.mastertable} --viz_dir_full {params.viz_dir_full} --viz_dir_relative {params.viz_dir_relative} --output_html {output.html_file}
        '''

rule casR_clustering:
    '''
    Creating clusters of everything annotated as CasR by CCTyper.
    Here we concatenate all CasRs, so we are not using loci-based wildcards.
    Instead we use output from concatenate_locus_viz as an indirect signal that all CasRs have been generated
    '''
    input:
        viz_done = rules.concatenate_locus_viz.output,
    output:
        clustered_CasR = base_path + "/casR_clustering/casR_clustered.faa",
        concatenated_casR = base_path + "/casR_clustering/casR_concatenated.faa",
    conda:
        "envs/trees.yaml"
    params:
        casR_wildcarded = base_path + "/14_cctyper_gene_locations/*/*casR.faa",
        cutoff = 0.4,
        n = 2,
    threads: thread_ultrasmall
    shell:
        '''
        cat {params.casR_wildcarded} > {output.concatenated_casR}
        cd-hit -i {output.concatenated_casR} -o {output.clustered_CasR} -c {params.cutoff} -n {params.n} -d 0 -M 16000 -T {threads}
        '''

#### Phage genomes ####

def aggregate_millard_phage_RN_analysis(wildcards):
    '''
    This function is used to aggregate the outputs of the millard_phage_RN_analysis rule.
    If not working, create a checkpoint to generate wildcards first.
    '''
    checkpoint_output = checkpoints.divide_millard_phages_to_folders.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{phage}/done.txt")).phage
    return expand(base_path + "/P3_hmm_analysis/{phage}/{phage}_cATyper_results.tsv", phage=cvals)


checkpoint divide_millard_phages_to_folders:
    '''
    Takes in the genome Millard genomes.fa file and creates a folder and fasta for each phage in the fasta.
    The headers are in the format >phage_name description. We only want to keep the phage_name for the folder.
    '''
    output: directory(base_path + "/P1_genomes")
    input:
        millard_fa = millard_fa, #nt multifasta of all phage genomes
        millard_proteins = millard_proteins #aa multifasta of all phage genones. Headers are >phage_name_runningnumber
    params:
        outdir = base_path + "/P1_genomes",
    threads: thread_ultrasmall
    run:
        from Bio import SeqIO
        import os
        print("Starting divide_millard_phages_to_folders")
        
        # Create a dictionary to store sequences for each phage
        phage_sequences = {}
        phage_protein_sequences = {}
        
        # Parse the millard_fa file and store sequences in the dictionary
        for record in SeqIO.parse(input.millard_fa, "fasta"):
            phage_name = record.id.split(" ")[0]
            if phage_name not in phage_sequences:
                phage_sequences[phage_name] = []
            phage_sequences[phage_name].append(record)
        
        # Parse the millard_proteins file and store sequences in the dictionary
        for record in SeqIO.parse(input.millard_proteins, "fasta"):
            phage_name = record.id.split("_")[0]
            if phage_name not in phage_protein_sequences:
                phage_protein_sequences[phage_name] = []
            phage_protein_sequences[phage_name].append(record)
        
        # Write sequences to files for each phage
        for phage_name, sequences in phage_sequences.items():
            phage_dir = f"{params.outdir}/{phage_name}"
            os.makedirs(phage_dir, exist_ok=True)
            
            # Write nucleotide sequences to file
            with open(f"{phage_dir}/{phage_name}.fna", "w") as out:
                SeqIO.write([seq for seq in sequences if seq.id.split(" ")[0] == phage_name], out, "fasta")
            
        for phage_name, sequences in phage_protein_sequences.items():
            phage_dir = f"{params.outdir}/{phage_name}"
            os.makedirs(phage_dir, exist_ok=True)

            # Write protein sequences to file
            with open(f"{phage_dir}/{phage_name}_proteins.faa", "w") as out:
                SeqIO.write([seq for seq in sequences if seq.id.split("_")[0] == phage_name], out, "fasta")
            
            # Write a "done" file in every phage folder
            with open(f"{phage_dir}/done.txt", "w") as out:
                out.write("done")

rule millard_phage_RN:
    '''
    Searches for ring nucleases in the Millard phage proteomes using HMMER and
    our pre-made RN hmm profiles similar to rule ring_nucleases
    '''
    input:
        phage_proteins = base_path + "/P1_genomes/{phage}/{phage}_proteins.faa"
    output:
        temp_rows = base_path + "/P2_hmm/{phage}/{phage}_RN_temp.tsv",
        hmm_rows = base_path + "/P2_hmm/{phage}/{phage}_RN_hmm.tsv"
    params:
        rn_db = ring_nuclease_db,
        evalue = "1e-8",
        temp_hmm = base_path + "/P2_hmm/{phage}/{phage}_RN_hmm.temp"
    conda: "envs/hmmer.yaml"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P2_hmm/{phage}/logs/RN_search.out",
        err = base_path + "/P2_hmm/{phage}/logs/RN_search.err"
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.rn_db} {input.phage_proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.phage}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.phage}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

rule millard_phage_RN_analysis:
    '''
    Analyses the ring nuclease hits in the Millard phage proteomes by reusing some code from rule ring_nucleases_analysis
    and phage-specific scripts too.
    '''
    input:
        hmm_rows = rules.millard_phage_RN.output.hmm_rows,
        proteins = base_path + "/P1_genomes/{phage}/{phage}_proteins.faa"
    output:
        catyper = base_path + "/P3_hmm_analysis/{phage}/{phage}_cATyper_results.tsv",
        hmm_targets = base_path + "/P3_hmm_analysis/{phage}/{phage}_RN_hmm_targets.tsv",
        effector_to_protein = base_path + "/P3_hmm_analysis/{phage}/{phage}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/P3_hmm_analysis/{phage}/{phage}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/P3_hmm_analysis/{phage}/{phage}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/P3_hmm_analysis/{phage}",
        hmm_msa_folder = ring_nuclease_folder + "/profiles"
    conda: "envs/hmmer.yaml"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P3_hmm_analysis/{phage}/logs/RN_analysis.out",
        err = base_path + "/P3_hmm_analysis/{phage}/logs/RN_analysis.err"
    shell:
        '''
        echo "Running ring nuclease analysis for phages" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/phage_RN_analyzer.py --sample {wildcards.phage} --output_folder {params.outdir} --proteins_fasta {input.proteins} --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "ring_nucleases" --effector_plot_data {output.plottable_effector_positions} --ring_nuclease True 2> {log.err} 1> {log.out}
        '''


rule concatenate_millard_phage_RN_analysis:
    input:
        files = aggregate_millard_phage_RN_analysis
    output: base_path + "/P3_hmm_analysis/phages_RN_hmm.tsv"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P3_hmm_analysis/logs/phages_RN_hmm.out",
        err = base_path + "/P3_hmm_analysis/logs/phages_RN_hmm.err"
    params:
        temp_file = base_path + "/P3_hmm_analysis/phages_RN_hmm.temp"
    shell:
        """
        # Get the header from the first file
        header=$(head -n 1 {input.files[0]})

        # Write the header to the output file
        echo "$header" > {output}

        # Concatenate files without header and append to the output file
        for file in {input.files}; do
            tail -n +2 "$file" >> {output}
        done
        """

rule phage_cas10:
    '''
    Searches the phage proteomes for Cas10s using the Cas10 hmm profiles
    '''
    input:
        phage_proteins = base_path + "/P1_genomes/{phage}/{phage}_proteins.faa"
    output:
        temp_rows = base_path + "/P4_cas10hmm/{phage}/{phage}_Cas10_temp.tsv",
        hmm_rows = base_path + "/P4_cas10hmm/{phage}/{phage}_Cas10_hmm.tsv"
    params:
        cas10_db = cas10_db,
        evalue = "1e-8",
        temp_hmm = base_path + "/P4_cas10hmm/{phage}/{phage}_Cas10_hmm.temp"
    conda: "envs/hmmer.yaml"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P4_cas10hmm/{phage}/logs/Cas10_search.out",
        err = base_path + "/P4_cas10hmm/{phage}/logs/Cas10_search.err"
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.cas10_db} {input.phage_proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.phage}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.phage}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

def aggregate_phage_cas10(wildcards):
    checkpoint_output = checkpoints.divide_millard_phages_to_folders.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{phage}/done.txt")).phage
    return expand(base_path + "/P4_cas10hmm/{phage}/{phage}_Cas10_hmm.tsv", phage=cvals)

rule phage_cas10_concatenate:
    input: aggregate_phage_cas10
    output: base_path + "/P4_cas10hmm/phages_Cas10_hmm.tsv"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P4_cas10hmm/logs/phages_Cas10_hmm.out",
        err = base_path + "/P4_cas10hmm/logs/phages_Cas10_hmm.err"
    params:
        temp_file = base_path + "/P4_cas10hmm/phages_Cas10_hmm.temp"
    shell:
        """
        # Get the header from the first file
        header=$(head -n 1 {input.files[0]})
        
        # Concatenate files without header
        tail -q -n +2 {input.files} > {params.temp_file}

        # Add header to the concatenated file
        echo "$header" > {output}

        # Concatenate the rest of the files
        cat {params.temp_file} >> {output}

        # Remove the temporary file
        rm {params.temp_file}
        """

#finish writing analysis if any cas10s are found
# rule phage_cas10_analysis:
#     '''
#     Takes the hmm hits from rule phage_cas10 and creates simplified tables similar to
#     ring nuclease rules before
#     '''
#     input:
#         hmm_rows = rules.phage_cas10.output.hmm_rows,
#         proteins = base_path + "/P1_genomes/{phage}/{phage}_proteins.faa"
#     output:
#         catyper = base_path + "/P5_cas10_analysis/{phage}/{phage}_cATyper_results.tsv",
#         hmm_targets = base_path + "/P5_cas10_analysis/{phage}/{phage}_Cas10_hmm_targets.tsv",
#         effector_to_protein = base_path + "/P5_cas10_analysis/{phage}/{phage}_effector_to_protein.tsv",
#         protein_to_effector = base_path + "/P5_cas10_analysis/{phage}/{phage}_protein_to_effector.tsv",
#         plottable_effector_positions = base_path + "/P5_cas10_analysis/{phage}/{phage}_plottable_effector_positions.tsv",
#     params:
#         outdir = base_path + "/P5_cas10_analysis/{phage}",
#         hmm_msa_folder = cas10_folder + "/profiles"
#     conda: "envs/hmmer.yaml"
#     threads: thread_ultrasmall
#     log:
#         out = base_path + "/P5_cas10_analysis/{phage}/logs/Cas10_analysis.out",
#         err = base_path + "/P5_cas10_analysis/{phage}/logs/Cas10_analysis.err"
#     shell:
#         '''
#         echo "Running Cas10 analysis for phages" >> {log.out}
#         echo "Listing targets" >> {log.out}
#         ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
#         python scripts/phage_RN_analyzer.py --sample {wildcards.phage} --output_folder {params.outdir} --proteins_fasta {input.proteins} --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "Cas10" --effector_plot_data {output.plottable_effector_positions} --ring_nuclease False 2> {log.err} 1> {log.out}
#         '''


rule spacers:
    input:
        done = rules.CRISPRCasTyper.output,
        renaming_done = rules.CRISPRCasTyper_rename.output
    params:
        fastafolder = base_path + "/07_cctyper/{j}/spacers",
        fastafile = base_path + "/07_cctyper/{j}/spacers/",
        crisprresultfolder = base_path + "/07_cctyper/{j}/"
    output:
        base_path + "/03_spacers/all/{j}_spacers.csv"
    run:
        from Bio import SeqIO
        import os.path
        from os import path
        import csv

        crispr_loci = [] #create a crispr_loci list. This is just a list of cctyper loci.
        crispr_results = str(params.crisprresultfolder + "crisprs_all.tab")
        print("Opening " + crispr_results)
        if os.path.isfile(crispr_results):
            with open(crispr_results) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                line_count = 0
                for row in csv_reader:
                    if line_count == 0:
                        #print(f'Column names are {", ".join(row)}')
                        line_count += 1
                    else:
                        #print(", ".join(row))
                        #print(row[11])
                        if row[11] == "True":
                            #print(row[1] + " trusted")
                            crispr_loci.append(row[1])
                        line_count += 1
                #print(f'Processed {line_count} lines.')
                #print("Trusted loci: " +  str(crispr_loci))
        else:
            print("No CRISPR-Cas loci found")

        fastafolder_object = Path(params.fastafolder) #store path to spacer folder as object
        shell("touch '{base_path}/03_spacers/all/{wildcards.j}_spacers.csv'") #create spacer .csv file for every sample
        if len(crispr_loci) > 0: #check if there are any loci in the trusted list
            print("------------- Spacers found in  " + str(wildcards.j) + "!")
            spacers_dict = {}
            for filename in os.listdir(fastafolder_object): #iterate through all fastas
                filename_reduced = re.search('(.*)\.[^.]+$', filename).group(1)
                #print("Reduced locus name: " + str(filename_reduced))
                #print(crispr_loci)
                #print(filename)
                if filename_reduced in crispr_loci:
                    #print("Match for file " + str(filename) + " in trusted loci")
                    for seq_record in SeqIO.parse(str(fastafolder_object) + '/' + filename, 'fasta'): #iterate through all entries in fasta
                        spacers_dict[seq_record.id] = seq_record.seq #add entry as dictionary entry
            with open(str(output), 'w') as f:
                for key in spacers_dict.keys():
                    f.write('%s,%s\n'%(key + "@" + str(wildcards.j),spacers_dict[key]))
        else:
            print("No spacers found in " + str(wildcards.j))

def aggregate_spacers(wildcards):
    '''
    Returns all spacers from all samples
    '''
    checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
    print(checkpoint_output)
    ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    print("Aggregating spacers from trusted locus observations")
    return expand(base_path + "/03_spacers/all/{j}_spacers.csv", j=ivals)

rule combined_spacers:
    '''
    This rule aggregates all spacer.csv files into one big .csv collection. The expand function is used for
    the concatenation, and it uses the list returned by the aggregate_IDs function.
    '''
    input: aggregate_spacers
    output: base_path + "/03_spacers/combined_spacers.csv"
    shell:
        """
        cat {input} > {output}
        """

rule spacers_blast_convert_spacers:
    """
    This converts the concatenated .csv into fasta format for blasting
    """
    input:
        rules.combined_spacers.output #the concatenated csv file
    output:
        fasta = base_path + "/03_spacers/combined_spacers.fasta",
        formatted_spacers = base_path + "/03_spacers/formatted_spacers.csv"
    threads: thread_ultrasmall
    shell:
        """
        sed 's/^/>/' {input} > {output.formatted_spacers}
        tr "," "\n" < {output.formatted_spacers} > {output.fasta}
        """

rule spacers_blast_runblast:
    '''
    Blasts spacers against the phage database.
    Also adds headers to blast result file.

    Remember to make a blast db of the phage genomes first.
    '''
    input:
        fasta = rules.spacers_blast_convert_spacers.output.fasta #the concatenated csv file
    output:
        blast = base_path + "/P5_spacer_to_phage_blast/blast.out",
    params:
        blast_db = millard_fa,
    conda: "envs/blast.yaml"
    threads: thread_hogger
    shell:
        '''
        blastn -num_threads {threads} -db {params.blast_db} -task blastn-short -query {input.fasta} -outfmt "6 qseqid qlen sseqid  stitle staxids saccver sstart send evalue length pident mismatch gapopen gaps sstrand" -evalue 0.01 > {output.blast}
        echo -e 'Query_id\tQuery_length\tSubject_id\tSubject_sci_title\tSubject_taxID\tSubject_accession_ID_version\tSubject_start\tSubject_end\tEvalue\tLength\tPercentage_identical\tMismatches\tOpen_gaps\tGaps\tSubject_strand' | cat - {output.blast} > blast.cat && mv blast.cat {output.blast}      
        '''

def aggregate_allHosts(wildcards):
    '''
    Returns all hosts after the filtering by cas10
    '''
    checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
    ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    return expand(base_path + "/03_postfiltering_genome_wildcards/{j}.txt", j=ivals)


#### RN trees and annotations for Crn4 paper ####
def aggregate_rn_wildcards(wildcards):
    checkpoint_output = checkpoints.rn_wildcarder.get(**wildcards).output[0]
    rn_vals = glob_wildcards(os.path.join(checkpoint_output,"{rn}.rn")).rn
    return expand(base_path + "/RN_01_wildcards/{rn}.rn", rn=rn_vals)

def aggregate_rn_trees(wildcards):
    checkpoint_output = checkpoints.rn_wildcarder.get(**wildcards).output[0]
    rn_vals = glob_wildcards(os.path.join(checkpoint_output,"{rn}.rn")).rn
    return expand(base_path + "/RN_04_tree/{rn}_tree.txt", rn=rn_vals)

checkpoint rn_wildcarder:
    '''
    This rule creates wildcards for each previosly known rn.
    Wildcards are used downstream when creating phylogenetic trees etc for each rn.
    '''
    output:
        directory(base_path + "/RN_01_wildcards")
    params:
        outdir = base_path + "/RN_01_wildcards",
    threads: thread_small
    run:
    #make the directory
        import os
        os.makedirs(params.outdir, exist_ok=True)
        rns = ["crn1", "crn2", "crn3", "crn4a", "crn4b", "csx15", "csx16", "csx20"]
        for rn in rns:
            print(rn)
            with open(f'{params.outdir}/{rn}.rn', 'w') as f:
                f.write(rn)


rule rn_fetcher:
    '''
    Fetches protein sequences for an effector from catyper_analysis outputs
    '''
    input:
        aggregate_rn_wildcards,
  #      rn_analysis_done = rules.concatenate_ring_nucleases_analysis.output, #this indicates we can safely probe the folders for rn sequences
    output:
        multifasta = base_path + "/RN_02_protseq/{rn}.faa",
    threads: thread_ultrasmall
    shell:
        '''
        echo "Fetching rn sequences for {wildcards.rn}"
        rn_name={wildcards.rn}
        touch {output.multifasta}
        find {base_path}/71_ring_nucleases_analysis/ -name \"*${{rn_name}}.faa\" -exec cat {{}} \; >> {output.multifasta}
        '''

rule rn_align:
    '''
    This rule copies logic from the effector_align rule
    '''
    input: rules.rn_fetcher.output.multifasta
    output: base_path + "/RN_03_aligment/{rn}.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/RN_03_aligment/logs/{rn}.out",
        err = base_path + "/RN_03_aligment/logs/{rn}.err"
    threads: thread_small
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out} 
        '''


rule rn_tree:
    input: rules.rn_align.output
    output: base_path + "/RN_04_tree/{rn}_tree.txt",
    conda: "envs/trees.yaml"
    threads: thread_small
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule concatenate_rn_wildcards:
    input:
        trees = aggregate_rn_trees,
        checkpoint = aggregate_rn_wildcards
    output: base_path + "/rns_finished.txt"
    threads: thread_small
    shell:
        '''
        cat {input} > {output}
        '''

#### Crn4 specific rules ####

rule enchanced_gff:
    '''
    Adds annotations from cctyper and other sources to the original genomic gff.
    Used when plotting locus visualisations of Crn4
    '''
    input:
        #genome_gff = rules.crispr_locus_proteins.output.crispr_locus_gff,
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        cctyper_plottable = rules.cctyper_gene_locations.output.cctyper_gene_locations_plottable,
        plottable_effector_positions = rules.validated_new_effectors_analysis.output.plottable_effector_positions,
        known_effectors = rules.cATyper_analysis.output.plottable_effector_positions,
        ring_nucleases = rules.ring_nucleases_analysis.output.plottable_effector_positions
    output:
        visualisation = base_path + "/gff_enhancer/{c}/{c}_enhanced.gff"
    conda:
        "envs/locus_visualiser.yaml"
    params:
        outdir = base_path + "/gff_enhancer/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        cctyper_folder = base_path + "/07_cctyper",
    log:
        out = base_path + "/gff_enhancer/logs/{c}.out",
        err = base_path + "/gff_enhancer/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/gff_enhancer.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --cctyper_folder {params.cctyper_folder} --validated_effectors {input.plottable_effector_positions} --cctyper_protein_table {input.cctyper_plottable} --known_effector_table {input.known_effectors} --ring_nucleases {input.ring_nucleases}  2> {log.err} 1> {log.out}
        '''

def aggregate_gff_enhancer(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/gff_enhancer/{c}/{c}_enhanced.gff", c=cvals)

rule concatenate_gff_enhancer:
    input: aggregate_gff_enhancer
    output: base_path + "/gff_enhancer/gff_enhancer.done"
    threads: thread_ultrasmall
    shell:
        '''
        touch {output}
        '''

def aggregate_crn4_genome_subsets(wildcards):
    checkpoint_output = checkpoints.create_crn4_genomes_wildcards.get(**wildcards).output[0]
    crn4_genome_vals = glob_wildcards(os.path.join(checkpoint_output,"{crn4_genome}/features.gff")).crn4_genome
    return expand(base_path + "/RN_06_crn4_gffs/{crn4_genome}/features_crn4_subset.gff", crn4_genome=crn4_genome_vals)

checkpoint create_crn4_genomes_wildcards:
    '''
    This checkpoint uses the .faa file of all crn4s to fetch the gff files of the genomes.
    Each crn4 protein sequence has a header that contains the genome accession number (e.g. >GCF_013177295.1_0__WP_235904890.1).
    The genome accession number is used to fetch the gff file of the genome (e.g. GCF_013177295.1.gff).
    '''
    input:
        crn4_faa = base_path + "/RN_02_protseq/crn4.faa"
    output:
        directory(base_path + "/RN_05_crn4_genomes")
    params:
        outdir = base_path + "/RN_05_crn4_genomes",
    threads: thread_small
    run:
        import os
        os.makedirs(params.outdir, exist_ok=True)
        with open(input.crn4_faa) as f:
            for line in f:
                if line.startswith(">"):
                    print("Processing line: " + line) #full line is e.g. GCF_013177295.1_1__WP_235904891.1
                    first_part = line.split("_")[0] #this is GCF or GCA
                    second_part = line.split("_")[1] #this is the accession number, e.g. 013177295.1
                    third_part = line.split("_")[2] #this is the locus number
                    locus = first_part + "_" + second_part + "_" + third_part # e.g. GCF_013177295.1_0
                    locus = locus.split(">")[1]
                    genome = locus.split("_")[0] + "_" + locus.split("_")[1] #e.g. GCF_013177295.1
                    #create the folder
                    os.makedirs(f"{params.outdir}/{locus}", exist_ok=True)
                    #original_gff_path = f"{base_path}/06_host_genomes/{genome}/{genome}_features.gff"
                    original_gff_path = f"{base_path}/gff_enhancer/{locus}/{locus}_enhanced.gff"
                    os.system(f"cp {original_gff_path} {params.outdir}/{locus}/features.gff")

                    protein_id = line.split("__")[1] #this is the protein id, e.g. WP_235904890
                    #ensure protein_id does not contain \n
                    protein_id = protein_id.strip()
                    print("Protein id: " + protein_id)
                    #make a text file with protein id in the genome folder
                    with open(f"{params.outdir}/{locus}/{protein_id}.txt", "w") as f:
                        f.write(protein_id)

rule subset_gffs_around_crn4:
    '''
    This rule subsets the gff files of the genomes around the crn4 proteins.
    '''
    input:
        gff = base_path + "/RN_05_crn4_genomes/{crn4_genome}/features.gff",
    output:
        subset_gff = base_path + "/RN_06_crn4_gffs/{crn4_genome}/features_crn4_subset.gff"
    params:
        gff_folder = base_path + "/RN_05_crn4_genomes/{crn4_genome}",
        subset_range = 20000
    log:
        out = base_path + "/RN_06_crn4_gffs/logs/{crn4_genome}.out",
    run:
        import os
        import gffutils
        db = gffutils.create_db(input.gff, dbfn=f"{params.gff_folder}/gff.db", force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        # list all .txt files in the folder
        txt_files = [f for f in os.listdir(params.gff_folder) if f.endswith('.txt')]
        protein_file = txt_files[0]
        with open(f"{params.gff_folder}/{protein_file}") as f:
            protein_id = f.readline().strip()

        #get the protein id from the protein id file
        print("Subsetting around " + str(protein_id), flush=True)
        try:
            protein_in_gff = db[protein_id]
        except:
            print("Did not find hit with protein name only. Adding cds- to the protein name.", flush=True)
            protein_id_csv = "cds-" + protein_id
            protein_in_gff = db[protein_id_csv]
        #defining surrounding region
        start = max(protein_in_gff.start - params.subset_range, 1)
        end = protein_in_gff.end + params.subset_range

        print("Surrounding region: " + str(start) + " - " + str(end), flush=True)

        subsetted_gff = db.region(region=(protein_in_gff.seqid, start, end), completely_within=False)

        with open (output.subset_gff, "w") as f:
            for feature in subsetted_gff:
                f.write(str(feature) + "\n")
            

rule concatenate_crn4_subsets:
    input: aggregate_crn4_genome_subsets
    output: base_path + "/crn4_genomes_subsetted"
    shell:
        '''
        cat {input} > {output}
        '''


rule final:
    input:
        tax = base_path + "/06_host_genomes/taxInfo.txt",
        tree_Cas10 = rules.Cas10_tree.output,
        type_iii_info = rules.combine_GGDD_HMM_to_mastertable.output.final_info_table,
        mastercombiner_final_info_table = rules.mastercombiner.output.final_info_table,
        concat_taxInfo = rules.concat_taxInfo.output,
        catyper_hmm = rules.concatenate_cATyper_hmm.output,
        catyper_analysis = rules.concatenate_cATyper_analysis.output,
        cas10_HD_faa = rules.Cas10_HD_hmm_maker.output.faa,
        cas10_hd_hmm = rules.Cas10_HD_hmmer.output,
        cas10_HD_hmm_merged = rules.merge_HD_hmm_with_info.output.merged_table,
        cas10_GGDD_faa = rules.Cas10_GGDD_hmm_maker.output.faa,
        cas10_GGDD_hmm = rules.Cas10_GGDD_hmmer.output,
        cas10_GGDD_hmm_merged = rules.merge_GGDD_hmm_with_info.output.merged_table,
        effector_commonness = rules.effector_commonness.output.effector_commonness_tsv,
        effector_scores_summary = rules.analyse_cATyper_effector_scores.output.effector_scores_summary,
        validated_effectors_scores_summary = rules.analyse_validated_new_effectors_scores.output.effector_scores_summary,
        concatenate_effector_wilcards = rules.concatenate_effector_wildcards.output,
        concatenate_cATyper_hmm_hhsuite = rules.concatenate_cATyper_hmm_hhsuite.output,
        parse_hhsuite = rules.concatenate_cATyper_hmm_hhsuite.output,
        parse_hhsuite_cogs = rules.concatenate_cATyper_hmm_hhsuite_cogs.output,
        concatenate_validated_new_effectors_analysis = rules.concatenate_validated_new_effectors_analysis.output,
        heatmap_known_validated_effectors = rules.heatmap_known_validated_effectors.output.effector_scores_summary,
        concatenate_locus_viz = rules.concatenate_locus_viz.output,
        html = rules.create_html_file.output.html_file,
        casR_cluster = rules.casR_clustering.output.clustered_CasR,
        ring_fusions = rules.ring_fusions.output.ring_fusions,
        ring_fusions_cas10 = rules.ring_nuclease_cas10_fusions.output.hmm_rows,
        validated_effectors_cas10_fusions = rules.validated_effectors_cas10_fusions.output.hmm_rows,
        known_effectors_cas10_fusions = rules.known_effectors_cas10_fusions.output.hmm_rows,
        ring_nucleases_fusions = rules.concatenate_ring_nuclease_fusions.output,
        concatenate_millard_phage_RN_analysis = rules.concatenate_millard_phage_RN_analysis.output,
        hostTable = rules.hostTable.output,
        concatenate_rn_wildcards = rules.concatenate_rn_wildcards.output,
        concatenate_crn4_subsets = rules.concatenate_crn4_subsets.output,
        concatenate_gff_enhancer = rules.concatenate_gff_enhancer.output,
        concatenate_ring_nucleases_analysis_no_length_lim = rules.concatenate_ring_nucleases_analysis_no_length_lim.output,
        concatenate_ring_nucleases_scores_no_length_lim = rules.concatenate_ring_nucleases_scores_no_length_lim.output,
    output: base_path + "/done"
    threads: thread_small
    shell:
        '''
        touch {output}
        mkdir -p {base_path}/R
        mkdir -p {base_path}/R/effector_trees

        mkdir -p {base_path}/R/effector_hmm
        mkdir -p {base_path}/R/effector_hmm/pfam
        mkdir -p {base_path}/R/effector_hmm/pdb
        mkdir -p {base_path}/R/effector_hmm/cog

        mkdir -p {base_path}/R/group4
        mkdir -p {base_path}/R/group4/pfam
        mkdir -p {base_path}/R/group4/pdb
        mkdir -p {base_path}/R/group4/cog
        mkdir -p {base_path}/R/group4/carfsaved
        
        cp {base_path}/12_Cas10_tree/cas10_tree.txt {base_path}/R
        cp {base_path}/06_host_genomes/taxInfo.txt {base_path}/R
        cp {base_path}/09_crispr_iii_CorA/loci/type_iii_info.tsv {base_path}/R
        
        cp {base_path}/13_effector_commonness/effector_commonness_master.tsv {base_path}/R
        cp {base_path}/mastertable_v2.tsv {base_path}/R

        cp -r {base_path}/45_effector_tree/* {base_path}/R/effector_trees

        cp -r {base_path}/43_effector_hmmer_analysis/*/*_sorted_filtered_mapped.tsv {base_path}/R/effector_hmm/pfam
        cp -r {base_path}/42_effector_hhsuite/*/*_parsed.tsv {base_path}/R/effector_hmm/pdb
        cp -r {base_path}/42_effector_hhsuite_cogs/*/*_parsed_cogs.tsv {base_path}/R/effector_hmm/cog

        '''