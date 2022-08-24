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
from opyrnd.jobs import JobRun, JobType
from opyrnd.utils import md5_from_file
from opyrnd.exceptions import OpyrndException,DoesNotExist,CouldNotConnect
"""
Deletes a All Entities and WorkFiles associated with a JobRun
"""



JOB_RUN_LABEL = 21
cfg = {'api_base_url': 'http://operend.bu.edu/api',"api_token_secret":'SOMETOKENSECRET','verify_https':'false'}
ApiBacked.configure_from_dict(cfg)
jr = JobRun.get_by_system_id(JOB_RUN_LABEL)
out_ids = []
for st in jr.subtasks:

    subtask_outfile_groups = st['outputWorkFileIds']
    print(subtask_outfile_groups)
    for og in subtask_outfile_groups:
        st_outfiles = subtask_outfile_groups[og]

        out_ids.extend(st_outfiles)
        print(len(st_outfiles))
print(len(out_ids))
print("Workfile Deletes")
deleted_wfs = 0
if hasattr(jr,'allGeneratedWorkFileIds'):
    print(len(jr.allGeneratedWorkFileIds))
    for wfid in jr.allGeneratedWorkFileIds:
        try:
            wf = WorkFile.get_by_system_id(wfid)
            d = wf.delete()
            deleted_wfs = deleted_wfs + 1
            print(f"DELETED {wfid} ")
        except DoesNotExist as de:
            print(f"{wfid} does not exist")
        except OpyrndException as oe:
            print(oe)
else:
    print("No WorkFiles to Delete in allGeneratedWorkFileIds")  
    sys.exit()

#all of the dangling entities
print("Entity Deletes")
deleted_ents = 0
for wfid in jr.allGeneratedWorkFileIds:
    print(wfid)
    e = Entity.get_by({'result_file': str(wfid)})
    print(e)
    if len(e) == 1:
        eid = e[0].entity_id
        d = Entity.delete_by_id(eid)
        deleted_ents = deleted_ents + 1
        print(f"DELETED {eid} {d}")
print("DONE")
print(f"WFS {deleted_wfs}")
print(f"ENTS {deleted_ents}")
