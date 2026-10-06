"""从已锁定的分子来源导出浏览器数据；不读取或导出患者表。"""
from pathlib import Path
import sys,json,hashlib,gzip

HERE=Path(__file__).resolve().parents[1]
CORE=HERE.parent/'interpreter'
sys.path.insert(0,str(CORE/'src'))
from tp53_interpreter.engine import Interpreter,clean

def main():
    e=Interpreter(); data=HERE/'docs/data'; data.mkdir(parents=True,exist_ok=True)
    fields=['kato_class','consensus_class','giacomelli_lof','dne_lof_class_nci','kotler_class_curated',
            'giacomelli_class_curated','clingen_preliminary_functional_code','alphamissense_score','alphamissense_class']
    variants={}; structures={}; specs=None; sources={}
    for p in e.quantifier.am.index:
        r=e.interpret(p); f=r['function']; axes=r['profile']['axes']
        if specs is None:
            specs=[{k:a[k] for k in ['metric','label','unit','reference_n']} for a in axes]
            # 首个变异不一定被全部实验覆盖；参考分母来自整个锁定来源，不能借首行的缺失。
            units=['kcal/mol','log10(mutant Kd / WT Kd)','percent of WT','Z-score','Z-score','Z-score','relative fitness score']
            for i,spec in enumerate(specs):
                spec['unit']=units[i]
                spec['reference_n']=[23,20][i] if i<2 else len(e.quantifier.populations[spec['metric']])
        variants[p[2:]]={'raw':[a['raw_value'] if a['raw_value'] is not None else a['raw_range'] for a in axes],
            'ranks':[a['percentile'] if a['percentile'] is not None else a['percentile_range'] for a in axes],
            'promoters':list(f['promoters_percent_wt'].values()),'classes':[f[k] for k in fields],
            'conflicts':r['conflicts'],
            'kotler_records':[[rid,v] for card in r['quantification']['source_assessments'] for m in card.get('measurements',[])
                               if m['metric']=='kotler_rfs' for rid,v in zip(m['source_record_ids'],m['record_values'])]}
        s=r['structure']
        if s['categories'] or s['residue_roles'] or s['measurements']:structures[p[2:]]=s
        for card in r['quantification']['source_assessments']:
            sources[card['evidence_family']]={k:card.get(k) for k in ['assay_context','assessment_kind','source']}
        if len(variants)%1000==0:print('exported',len(variants),flush=True)
    mutations=e.tables['mutation']; annotations={}
    cols=['SpliceAI_DS_AL','SpliceAI_DS_AG','SpliceAI_DS_DG','SpliceAI_DS_DL']
    for c,frame in mutations.groupby('c_description'):
        if not isinstance(c,str) or not c.startswith('c.'):continue
        values=frame[[v for v in cols if v in frame]].stack().dropna().tolist()
        annotations[c]={'protein':sorted(set(frame.ProtDescription.dropna())),
                        'spliceai_max':float(max(values)) if values else None,
                        'spliceai_complete':bool(set(cols)<=set(frame.columns) and frame[cols].notna().all().all())}
    ref=e.ref
    reference={'sequence':e.seq,'cds':ref['cds_v6'],
        'tails':{str(v):ref[f'transcript_sequence_v{v}'][ref[f'cds_start_transcript_v{v}']-1:] for v in [5,6]},
        'exons':{str(v):ref[f'exons_v{v}'] for v in [5,6]}}
    package={'version':'0.4.0','snapshot_id':e.snapshot_id,'date':'2026-10-06','scope':'public molecular sources only',
             'class_fields':fields,'axis_specs':specs,'reference':reference,'sources':sources,
             'variants':variants,'structures':structures,'coding_annotations':annotations}
    path=data/'public.json';path.write_text(json.dumps(clean(package),ensure_ascii=False,separators=(',',':'),allow_nan=False))
    metadata={'version':'0.4.0','molecular_snapshot':e.snapshot_id,'source_versions':e.versions,
              'variant_n':len(variants),'structural_measured_n':sum(bool(v['measurements']) for v in structures.values()),
              'public_json_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'public_json_bytes':path.stat().st_size,
              'exported_fields':['public experimental ranks and values','public curated classes','public sequence and exon boundaries'],
              'patient_data_exported':False}
    (data/'manifest.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2))
    (data/'public.json.gz').write_bytes(gzip.compress(path.read_bytes(),mtime=0))
    fixtures=HERE/'tests/fixtures';fixtures.mkdir(parents=True,exist_ok=True)
    for name in ['all_coding_snv_v03.csv','all_single_base_indels_v03.csv']:
        source=CORE/'outputs/interpreter_v03/02_清洗数据'/name
        (fixtures/name).write_bytes(source.read_bytes())
    small=[]
    for example_path in (CORE/'outputs/interpreter_v03/05_分析报告/逐变异').glob('*.json'):
        r=json.loads(example_path.read_text());small.append({k:r[k] for k in ['query','normalization','profile','non_missense'] if k in r})
    (fixtures/'representative.json').write_text(json.dumps(small,ensure_ascii=False,separators=(',',':')))
    print(json.dumps({'variants':len(variants),'bytes':path.stat().st_size,'snapshot':e.snapshot_id},ensure_ascii=False))

if __name__=='__main__':main()
