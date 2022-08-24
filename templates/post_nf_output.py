#!/usr/bin/env python3
from genericpath import exists
import os 
import sys
import csv
import pprint
import json
import argparse
from opyrnd.api_backed import ApiBacked
from opyrnd.workfile import WorkFile
from opyrnd.entities import Entity
from opyrnd.entity_classes import EntityClass
from opyrnd.jobs import JobRun,JobType
from opyrnd.utils import md5_from_file
from datetime import datetime as dt

OUT_DIR = "/rprojectnb/pulmseq/daupipe/OUTPUT/"
"""
Copied from a command line version in
the hand_scripts dir
"""

INPUT_FILE = "$input_tsv_file"
OUT_DIR = "$output_dir"
JOB_TYPE_NAME = "$job_type_name"
SEQ_READ_SET = "$seq_read_set"
API_BASE_URL = "$opapi_base_url"
API_VERIFY_HTTPS = "$opapi_verify_https"
API_STORAGE_LOCATION = "$opapi_storage_location"
JOB_RUN_LABEL = "$job_run_label"
TEST_RUN = True

SAMPLE_ISR_DICT = {}
def fq_ents_from_config(config_file,sample_id):
    """
    Returns a list containing the IDs of all Entities that were found
    for the fastqs of a given sample -- YES ALL of the
    Entities as there may be more than one. This assumes that there
    are no dupes for the FastQs.
    
    @param: config_file the config file used in the nextflow run 
    @param: the id of the sample if only one set of fastqs is needed
    otherwise if not provided, all ents will be returned
    """
    if sample_id and sample_id in SAMPLE_ISR_DICT:
        return SAMPLE_ISR_DICT[sample_id]

    ent_ids = []
    fq_ids = []
    with open(config_file) as conf_file:
        conf_reader = csv.DictReader(conf_file, delimiter="\t")
        for row in conf_reader:
            row_sample = row['SAMPLE_ID'] 
            #print(row)
            #print(row['SAMPLE_ID'])
            #print(row['R1'])
            #print(row['R2'])
            #print(os.path.isfile(row['R1']))
            #print(os.path.isfile(row['R2']))
            if sample_id and not row_sample == sample_id:
                continue
            md51 = md5_from_file(row['R1'])
            md52 = md5_from_file(row['R2'])
            wf1 = WorkFile.get_by_hash(md51)
            wf2 = WorkFile.get_by_hash(md52)
            #just use the first one we get - hopefully no dupes!
            if len(wf1)>0 and len(wf2)>0 :
                print("THE WF IDs")
                print(str(wf1[0].systemId))
                print(str(wf2[0].systemId))
                fq_ids.append((wf1[0].systemId))
                fq_ids.append((wf2[0].systemId))
                # yes may really be more than one ent returned
                ents = Entity.get_by({'_class':'IlluminaSeqReads','FASTQ':str(wf1[0].systemId),'FASTQ':str(wf2[0].systemId)})
                if len(ents) > 0:
                    print(ents[0].to_dict())
                else:
                    print("NOTHING")
                print("MY ENTS")
                print(ents)   
                ent_ids.extend([e.entity_id for e in ents])
                #ent.__setitem__('illumina_seq_reads',ent[0]._entity_id)
                #ent.save()
    SAMPLE_ISR_DICT[sample_id] = (ent_ids,fq_ids)
    return ent_ids,fq_ids

def post_runSTAR1pass(input_rows):
    # a single file - {SAMPLE_ID}.1pass.SJ.out.tab
    """
    Post the file
    return {name="runSTAR1pass',comments="dunno",{outputWorkFileIds:[filed]}
    or .. just the file IDs - we can add the name and such when calling this 
    """
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    for row in input_rows:
        path = f"{row['INDIVIDUAL_ID']}/{row['SAMPLE_ID']}/Processing/Libraries/{row['LIBRARY_ID']}/{row['RG_ID']}/STAR_1Pass/{row['SAMPLE_ID']}.1pass.SJ.out.tab"
        path = os.path.join(OUT_DIR,path)
        file_paths.append(path)
        if not TEST_RUN:
            wf = WorkFile.post_from_file(path)
            print(f"Posted WorkFile {wf}")
            # The Entity uses lower case field names  
            ent_fields = dict((key.lower(), value) for (key, value) in row.items())
            ent_fields['next_flow_process'] = "STAR1pass"
            ent_fields["result_file"]=wf.systemId
            ent_fields["result_file_original_name"] = os.path.basename(path)
            ent_bodies.append(ent_fields)
            w_file_ids.append(wf.systemId)
            #print(os.stat(path))
    return w_file_ids,ent_bodies,file_paths

def post_runSTAR2pass(input_rows):
    # 3 files - {SAMPLE_ID}.bai,bam and bam.bai
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    for row in input_rows:
        #print("2 MY FILES")
        for ext in ["bam","bai","bam.bai"]:
            path = os.path.join(OUT_DIR,f"{row['INDIVIDUAL_ID']}/{row['SAMPLE_ID']}/{row['SAMPLE_ID']}.{ext}")
            #print(os.stat(path))
            file_paths.append(path)
            if not TEST_RUN:
                wf = WorkFile.post_from_file(path)
                
                print(f"Posted WorkFile {wf}")
                ent_fields = dict((key.lower(), value) for (key, value) in row.items())
                ent_fields['next_flow_process'] = "STAR2pass"
                ent_fields["result_file"]=wf.systemId
                ent_fields["result_file_original_name"] = os.path.basename(path)
                ent_bodies.append(ent_fields)
                w_file_ids.append(wf.systemId)
    return w_file_ids,ent_bodies,file_paths

def post_runRSEM(input_rows):
    # all files in the RSEM dir
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    #print("3 MY FILES")
    for row in input_rows:
        path = os.path.join(OUT_DIR,f"{row['INDIVIDUAL_ID']}/{row['SAMPLE_ID']}/RSEM")
        #print(path)
        for f in os.listdir(path):
            ent_fields = dict((key.lower(), value) for (key, value) in row.items())
            ent_fields['next_flow_process'] = "RSEM"
            file_paths.append(os.path.join(path,f))
            if not TEST_RUN:
                wf = WorkFile.post_from_file(os.path.join(path,f))
                print(f"Posted WorkFile {wf}")
                ent_fields["result_file"]=wf.systemId
                ent_fields["result_file_original_name"] = f
                ent_bodies.append(ent_fields)
                w_file_ids.append(wf.systemId)
                #print(f)
                #print(os.stat(os.path.join(path,f)))
    return w_file_ids,ent_bodies,file_paths

def post_runFastQC(input_rows):
    # a zip file for each end of a sample or - just everything in the dir
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    #print("3 MY FILES")
    for row in input_rows:
        path = f"{row['INDIVIDUAL_ID']}/{row['SAMPLE_ID']}/Processing/Libraries/{row['LIBRARY_ID']}/{row['RG_ID']}/FastQC"
        path = os.path.join(OUT_DIR,path)
        for f in os.listdir(path):
            file_paths.append(os.path.join(path,f))
            if not TEST_RUN:
                wf = WorkFile.post_from_file(os.path.join(path,f))
                print(f"Posted WorkFile {wf}")
                ent_fields = dict((key.lower(), value) for (key, value) in row.items())
                ent_fields['next_flow_process'] = "FastQC"
                ent_fields["result_file"]=wf.systemId
                ent_fields["result_file_original_name"] = f
                ent_bodies.append(ent_fields)
                w_file_ids.append(wf.systemId)
                #print(os.stat(os.path.join(path,f)))
    return w_file_ids,ent_bodies,file_paths

def post_runRSeQC(input_rows):
    # return {"name="runRSeQC",comments="dunno",{outputWorkFileIds:[wf_ids]}
    # lots of files in the dir
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    #print("4 MY FILES")
    for row in input_rows:
        path = f"{row['INDIVIDUAL_ID']}/{row['SAMPLE_ID']}/RSeQC" 
        path = os.path.join(OUT_DIR,path)
        for f in os.listdir(path):
            file_paths.append(os.path.join(path,f))
            if not TEST_RUN:
                wf = WorkFile.post_from_file(os.path.join(path,f))
                print(f"Posted WorkFile {wf}")
                ent_fields = dict((key.lower(), value) for (key, value) in row.items())
                ent_fields['next_flow_process'] = "RSeQC"
                ent_fields["result_file"]=wf.systemId
                ent_fields["result_file_original_name"] = f
                ent_bodies.append(ent_fields)
                w_file_ids.append(wf.systemId)
                #print(os.stat(os.path.join(path,f)))
    return w_file_ids,ent_bodies,file_paths

def post_runMultiQCFastq():
    # an html file and everything in the fastq_multiqc_data directory
    #print("5 MY FILES")
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    path = os.path.join(OUT_DIR,"Output/QC/Fastq")

    html_file = os.path.join(path,'fastq_multiqc.html')
    file_paths.append(html_file)
    if not TEST_RUN:
        wf = WorkFile.post_from_file(html_file)
        ent_fields = {"category":"QC",'next_flow_process':"MultiQCFastq","label":"Fastq","result_file":wf.systemId}
        ent_fields["result_file_original_name"] = 'fastq_multiqc.html'
        ent_bodies.append(ent_fields)
        #print(html_file)
        os.stat(html_file)
        w_file_ids.append(wf.systemId)
    data_path = os.path.join(path,'fastq_multiqc_data') 
    for f in os.listdir(data_path):
        #print(f)
        file_paths.append(os.path.join(data_path,f))
        if not TEST_RUN:
            wf = WorkFile.post_from_file(os.path.join(data_path,f))
            print(f"Posted WorkFile {wf}")
            ent_fields = {"category":"QC",'next_flow_process':"MultiQCFastq","label":"Fastq","result_file":wf.systemId}
            ent_fields["result_file_original_name"] =f
            ent_bodies.append(ent_fields)
            w_file_ids.append(wf.systemId)
            #print(os.stat(os.path.join(data_path,f)))
    return w_file_ids,ent_bodies,file_paths

def post_runMultiQCSample():
    # an html file and everything in the sample_multiqc_data directory
    #print("6 MY FILES")
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    path = os.path.join(OUT_DIR,"Output/QC/Sample")
    html_file = os.path.join(path,'sample_multiqc.html')
    file_paths.append(html_file)
    if not TEST_RUN:
        wf = WorkFile.post_from_file(html_file)
        
        print(f"Posted WorkFile {wf}")
        ent_fields = {"category":"QC",'next_flow_process':"MultiQCSample","label":"Sample","result_file":wf.systemId}
        ent_fields["result_file_original_name"] = 'sample_multiqc.html'
        ent_bodies.append(ent_fields)
        w_file_ids.append(wf.systemId)

        #print(html_file)
        os.stat(html_file)
    data_path = os.path.join(path,'sample_multiqc_data') 
    for f in os.listdir(data_path):
        #print(f)
        file_paths.append(os.path.join(data_path,f))
        if not TEST_RUN:
            wf = WorkFile.post_from_file(os.path.join(data_path,f))
            print(f"Posted WorkFile {wf}")
            ent_fields = {"category":"QC",'next_flow_process':"MultiQCSample","label":"Sample","result_file":wf.systemId}
            ent_fields["result_file_original_name"] =f
            ent_bodies.append(ent_fields)
            w_file_ids.append(wf.systemId)
            #print(os.stat(os.path.join(data_path,f)))
    return w_file_ids,ent_bodies,file_paths

def post_runSTARgenomeGenerate():
    # everything in the STAR/genomeGenerate dir
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    path=os.path.join(OUT_DIR,"Output/STAR/genomeGenerate")
    for f in os.listdir(path):
        file_paths.append(os.path.join(path,f))
        if not TEST_RUN:
            wf = WorkFile.post_from_file(os.path.join(path,f))
            print(f"Posted WorkFile {wf}")
            ent_fields = {"category":"STAR",'next_flow_process':"STARgenomeGenerate","label":"genomeGenerate","result_file":wf.systemId}
            ent_fields["result_file_original_name"] =f
            ent_bodies.append(ent_fields)
            w_file_ids.append(wf.systemId)
            #print(os.stat(os.path.join(path,f)))
    return w_file_ids,ent_bodies,file_paths

def post_createSE():
    # all the rds files - or just everything here
    path=os.path.join(OUT_DIR,"Output/Expression")
    print("THE PATH")
    print(path)
    w_file_ids = []
    ent_bodies = []
    file_paths = []
    for f in os.listdir(path):
        print(f)
        file_paths.append(os.path.join(path,f))
        if not TEST_RUN:
            wf = WorkFile.post_from_file(os.path.join(path,f))
            print(f"Posted WorkFile {wf}")
            w_file_ids.append(wf.systemId)
            ent_fields = {"category":"Expression",'next_flow_process':"createSE","label":"rds","result_file":wf.systemId}
            ent_fields["result_file_original_name"] =f
            ent_bodies.append(ent_fields)
            #print(os.stat(os.path.join(path,f)))
    return w_file_ids,ent_bodies,file_paths

def post_job_run():
    return 32

def print_file_paths(file_paths):
    #print(file_paths)
    for fp in file_paths:
        if os.path.exists(fp):
            print(f"{os.path.getsize(fp)}\t{fp}")
        else:
            print(f"-1\t{fp}")

def post_results(project_name,study_name,seq_read_set_name,input_tsv_ID):
    sub_tasks=[]
    sample_ent_bodies=[]
    coll_ent_bodies=[]
    file_paths = []
    # Post all Output Files
    # gather results of the createSE process
    cse = post_createSE()
    sub_tasks.append({"name":"createSE","comments":"one","outputWorkFileIds":{'outs':cse[0]}})
    coll_ent_bodies.append(cse[1])
    file_paths.extend(cse[2])
    # gather results of the MultiQCFastq process
    mqcf = post_runMultiQCFastq()
    sub_tasks.append({"name":"MultiQCFastq","comments":"one","outputWorkFileIds":{'outs':mqcf[0]}})
    coll_ent_bodies.append(mqcf[1])
    file_paths.extend(mqcf[2])
    # gather results of the MultiQCFastq process
    mqcs = post_runMultiQCSample()
    sub_tasks.append({"name":"MultiQCSample","comments":"one","outputWorkFileIds":{'outs':mqcs[0]}})
    coll_ent_bodies.append(mqcs[1])
    file_paths.extend(mqcs[2])
    # gather results of the STARgenomeGenerate process
    sgg = post_runSTARgenomeGenerate()
    sub_tasks.append({"name":"STARgenomeGenerate","comments":"one","outputWorkFileIds":{'outs':sgg[0]}})
    coll_ent_bodies.append(sgg[1])
    file_paths.extend(sgg[2])
    tsv_file = open(INPUT_FILE, "r")
    dict_reader = csv.DictReader(tsv_file,delimiter="\t")
    input_list = list(dict_reader)
    # gather results of the STAR1pass process
    sp1 = post_runSTAR1pass(input_list)
    sub_tasks.append({"name":"STAR1pass","comments":"one","outputWorkFileIds":{'outs':sp1[0]}})
    sample_ent_bodies.append(sp1[1])
    file_paths.extend(sp1[2])
    # gather results of the STAR2pass process
    sp2 = post_runSTAR2pass(input_list)
    sub_tasks.append({"name":"STAR2pass","comments":"one","outputWorkFileIds":{'outs':sp2[0]}})
    sample_ent_bodies.append(sp2[1])
    file_paths.extend(sp2[2])
    # gather results of the RSEM process
    rsm = post_runRSEM(input_list)
    sub_tasks.append({"name":"RSEM","comments":"one","outputWorkFileIds":{'outs':rsm[0]}})
    sample_ent_bodies.append(rsm[1])
    file_paths.extend(rsm[2])
    # gather results of the FastQC process
    fqc = post_runFastQC(input_list)
    sub_tasks.append({"name":"FastQC","comments":"one","outputWorkFileIds":{'outs':fqc[0]}})
    sample_ent_bodies.append(fqc[1])
    file_paths.extend(fqc[2])
    # gather results of the RSeQC process
    rqc = post_runRSeQC(input_list)
    sub_tasks.append({"name":"RSeQC","comments":"one","outputWorkFileIds":{'outs':rqc[0]}})
    sample_ent_bodies.append(rqc[1])
    file_paths.extend(rqc[2])
    #this is crap, need the real input workfiles some how
    # also need to set the output files
    # and need to set the job run ID or maybe datetime to disinguish between
    # different results - otherwise they'd all be the same
    ent_ids,fq_ids = fq_ents_from_config(INPUT_FILE,None)
    jr_dict = {'subtasks':sub_tasks,'jobTypeName':JOB_TYPE_NAME,'inputWorkFileIds':{'fq_files':fq_ids,'input_tsv':input_tsv_ID},'status':'COMPLETE'}
    #print(f"SAVING JR AS {json.dumps(jr_dict, indent = 4)}")
    #post the job run
    if TEST_RUN:
        print_file_paths(file_paths)
    if not TEST_RUN:
        jrun = JobRun(jr_dict)
        save_ouput = jrun.save()
        print("JUST SAVED JOB RUN")
        print(save_ouput)
        #print(f"SAVED JR {jrun}")
        #jr_json =json.dumps(jrun, indent = 4)
        #print(jr_json)
        # update ents with project and study
        #sample_ent_bodies is a list of lists
        #that contain dicts
        now = dt.now()
        for bods in sample_ent_bodies:
            for b in bods:
                #get the seq reads object for the sample
                ent_ids,fq_ids = fq_ents_from_config(INPUT_FILE,b['sample_id'])
                #gotta do something about the possibility of more than one Ent
                #for one of these records

                if len(ent_ids) > 0:
                    b['illumina_seq_reads'] = ent_ids[0]
                #ignore the run_date for now, the format could cause grief
                b.pop('run_date')
                b.pop('r1')
                b.pop('r2')
                b['seq_read_set_name'] = seq_read_set_name
                # no project or study now -- maybe delete this
                #b['project'] = project_name
                #b['study'] = study_name
                b['job_run_label'] = JOB_RUN_LABEL + "_" + now.strftime("%m_%y_%I:%M_%p")
                b['job_run'] = str(jrun.id)
                print("PSTING AN ENT")
                print(b)

                samp_ent = Entity(entity_class='RNASeqRes_PerSample',
                                values=b)
                print(samp_ent.to_dict())
                samp_ent.save()
                print(f"POSTED ENT {samp_ent}")
        for ebods in coll_ent_bodies:
            for e in ebods:
                #e['project'] = project_name
                #e['study'] = study_name
                b['seq_read_set_name'] = seq_read_set_name
                e['job_run_label'] = JOB_RUN_LABEL + "_" + now.strftime("%m_%y_%I:%M_%p")
                e['job_run'] = str(jrun.id)
                coll_ent = Entity(entity_class='RNASeqRes_PerJob',
                                values=e)
                coll_ent.save()
                print(f"POSTED ENT {coll_ent}")
        #Post a JobRun
        #Post the Ents
        #print(json.dumps(sample_ent_bodies, indent = 4))
        #print(json.dumps(coll_ent_bodies, indent = 4)) 
"""
for indivID,sampleIDs in samples.items():

    print(indivID)
    post_runSTAR1pass(indivID,sampleIDs)
    post_runSTAR2pass(indivID,sampleIDs)
    post_runRSEM(indivID,sampleIDs)
    post_runFastQC(indivID,sampleIDs)
    post_runRSeQC(indivID,sampleIDs)
    
"""
def pre_test():
    """
    Return a simple list of strings
    with error codes or None if all
    is well
    """
    errs = []
    # check that the job type exists
    if(not JobType.get_by_name(JOB_TYPE_NAME)):
        errs.append(f'Job Type {JOB_TYPE_NAME} Not Found in Core System')
    if(not EntityClass.get_by_name("RNASeqRes_PerSample")):
        errs.append('Entity Class RNASeqRes_PerSample Not Found in Core System')
    if (not EntityClass.get_by_name("RNASeqRes_Collective")):
        errs.append(f'Entity Class RNASeqRes_Collective Not Found in Core System')

    return errs

def run_as_template():
    """
    Run the output loader as a template from within nextflow. Expects
    the DOLLAR_SIGNvariables to be set in the nextflow process. The
    api_token_secret is set via a nextflow secret and must be retrieved
    via a env var: OPAPI_TOKEN_SECRET
    """

    global TEST_RUN
    global JOB_RUN_LABEL

    cfg = {'api_base_url': API_BASE_URL,
            'api_token_secret':os.getenv("OPAPI_TOKEN_SECRET"),
            'verify_https':API_VERIFY_HTTPS,
            'storage_location':API_STORAGE_LOCATION}
    
    ApiBacked.configure_from_dict(cfg)
    # no TEST_RUN within nextflow!
    TEST_RUN = False   
    input_wf_id = -999
    input_wf_id = WorkFile.post_from_file(INPUT_FILE).systemId

    #not using project or study at the moment
    post_results(None,None,SEQ_READ_SET,input_wf_id)

def run_as_standalone():
    """
    Run as standalone - well duh
    """
    global OUT_DIR
    global JOB_TYPE_NAME
    global INPUT_FILE
    global TEST_RUN
    global JOB_RUN_LABEL
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--operend_config_file", type=str,
                        help="The configuration ini file containing Operend configuration parameters",
                        required=True)
    parser.add_argument("-r","--results_directory", type=str,
                        help="The directory containing the results or publish directory of a pipeline run. This is where the files are read from",
                        required=True)

    parser.add_argument("-j","--job_type_name", type=str,
                        help="The name of the job type. Must be a valid Job Type that the user has access to within the Operend System",
                        required=True)
    
    parser.add_argument("-l","--job_run_label", type=str,
                        help="A text label to identify this specific job run. The datetime will be appended to the label when job run completes and is loaded",
                        required=True)
    
    parser.add_argument("-i","--input_file", type=str,
                        help="the full path to the tsv file containing the columns needed to run and store the results of the pipeline",
                        required=True)
    
    parser.add_argument("-s","--seq_read_set_name", type=str,
                        help="The name of the IlluminaSeqReadSet used as the source of the fastq input files",
                        required=True)
    parser.add_argument("-d","--dryrun",action='store_true',
                        help="if set will do a dry run and not store anything")
    
    args = parser.parse_args()
    
    args_good = True
    if not os.path.isfile(args.operend_config_file):
        print(f"Operend config file not found: {args.operend_config_file}")
        args_good = False
    else:
        ApiBacked.configure_from_file(args.operend_config_file)
    # test for the res dir
    if not os.path.isdir(args.results_directory):
        print(f"results directory not found: {args.results_directory}")
        args_good = False
    else:
        OUT_DIR = args.results_directory
    # test for job type name - better exist or we're SOL
    # this thing is basically hard coded to it 
    if JobType.get_by_name(args.job_type_name) is None:
        print(f"job type with specified name does not exist: {args.job_type_name}")
        args_good = False
    else:
        JOB_TYPE_NAME = args.job_type_name

    if not os.path.isfile(args.input_file):
        print(f"input_file not found: {args.input_file}")
        args_good = False
    else:
        INPUT_FILE = args.input_file
    
    

    if not Entity.get_by({'_class':'IlluminaSeqReadsSet','name':args.seq_read_set_name}):
        print(f"IlluminaSeqReadsSet not found: {args.seq_read_set_name}")
        args_good = False
        #hmm - why isnt seq read set name global?
    
    JOB_RUN_LABEL = args.job_run_label
    TEST_RUN = args.dryrun
    if not args_good:
        sys.exit()
    
    input_wf_id = -999
    if not TEST_RUN:
        input_wf_id = WorkFile.post_from_file(INPUT_FILE).systemId

    #not using project or study at the moment
    post_results(None,None,SEQ_READ_SET,input_wf_id)
    

if __name__ == "__main__":
    # if any of the template variables are set
    # to anything other than the template vals
    # we treat this as template file in a nextflow run
    # kinda hacky here -- should look for the dollar
    # sign but that would mean figuring out how to get
    # nextflow's template system to not replace it
    if (INPUT_FILE.find("input_tsv_file") == -1 or
        OUT_DIR.find("output_dir") == -1 or
        JOB_TYPE_NAME.find("job_type_name") == -1 or
        SEQ_READ_SET.find("seq_read_set") == -1 or
        API_BASE_URL.find("opapi_base_url") == -1 or
        API_VERIFY_HTTPS.find("opapi_verify_https") == -1 or
        API_STORAGE_LOCATION.find("opapi_storage_location") == -1
        ):
        run_as_template()
    else:
        run_as_standalone()

    sys.exit()

    