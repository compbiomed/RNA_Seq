#!/usr/bin/python3
from opyrnd.api_backed import ApiBacked
from opyrnd.entities import Entity
from opyrnd.workfile import WorkFile
cfg = {'api_base_url': 'http://operend.bu.edu/api',
	"api_token_secret":'SOME SECRET THING',
	'verify_https':'false'}
ApiBacked.configure_from_dict(cfg)
seq_read_sets = Entity.get_by({'_class':"IlluminaSeqReadsSet","name": "MCL2"})
#print(seq_read_sets[0].to_dict())
tsv_file = open("this_tsv_file_MCL2.tsv","w")
tsv_file.write('INDIVIDUAL_ID\tSAMPLE_ID\tLIBRARY_ID\tRG_ID\tPLATFORM_UNIT\tPLATFORM\tPLATFORM_MODEL\tRUN_DATE\tCENTER\tR1\tR2\n')
for sr_id in seq_read_sets[0].__getitem__('seqReads'):
	seq_read = Entity.get_by_id(sr_id)
	#print("SEQ READ")
	#print(seq_read.to_dict())
	individual_id=""
	sample_id=""
	library_id=""
	rg_id=""
	platform_unit=""
	platform=""
	platform_model=""
	run_date=""
	center=""
	r1=seq_read.__getitem__('FASTQ')[0]
	r2=seq_read.__getitem__('FASTQ')[1]
	wf1 = WorkFile.get_by_system_id(r1)
	wf2 = WorkFile.get_by_system_id(r2)
	r1_path = f"/rprojectnb/pulmseq/daupipe/PLINE_SCRATCH_SM/{wf1.originalName}"
	r2_path = f"/rprojectnb/pulmseq/daupipe/PLINE_SCRATCH_SM/{wf2.originalName}"
	# copy the workfiles
	#for wf in [wf1,wf2]:
	#		new_file = open(f"/rprojectnb/pulmseq/daupipe/PLINE_SCRATCH/{wf.originalName}", "wb")
	#		new_file.write(wf.open().read())
	for pu_id in seq_read.__getitem__('platformUnits'):
			pf_unit = Entity.get_by_id(pu_id)
			if len(platform_unit) > 0:
					platform_unit+=","
			platform_unit += pf_unit.__getitem__("flowcell")
			lanes = pf_unit.__getitem__('lanes')
			platform_unit +="."
			platform_unit += (''.join(str(l)+"," for l in lanes))[:-1]
			barcodes = pf_unit.__getitem__("barcodes")
			platform_unit +="." + '+'.join(barcodes)
			if len(rg_id) >0:
					rg_id += ","
			rg_id += pf_unit.__getitem__("flowcell") + "."+  ''.join(str(l)+"," for l in lanes)
			#get rid of that last comma!
			rg_id = rg_id[:-1]
			lib_id = pf_unit.__getitem__('library')
			seq_library = Entity.get_by_id(lib_id)
			#print(seq_library.to_dict())
			#print("SAMPLE")
			bs_id = seq_library.__getitem__('sample')
			bio_sample = Entity.get_by_id(bs_id)
			individual_id = sample_id = library_id = bio_sample.__getitem__('ID')
			#print(bio_sample.to_dict())
			run_id = pf_unit.__getitem__('run')
			seq_run = Entity.get_by_id(run_id)
			#print(seq_run.to_dict())
			inst_id = seq_run.__getitem__('instrument')
			instrument = Entity.get_by_id(inst_id)
			#print(instrument.to_dict())
			platform = instrument.__getitem__('manufacturer')
			platform_model = instrument.__getitem__('model')
			center = instrument.__getitem__('center')
			#print(instrument.to_dict())
	tsv_file.write(f"{individual_id}\t{sample_id}\t{library_id}\t{rg_id}\t{platform_unit}\t{platform}\t{platform_model}\t{run_date}\t{center}\t{r1_path}\t{r2_path}\n")
