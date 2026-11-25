"""Recompute historical descriptive panels from a hash-locked RNN table.

The manuscript curves retain row-level descriptive SEM; model-aware inference
is supplied by the separately included revised statistics module.
"""

import argparse,csv,json,os
from pathlib import Path
try:
    from .rnn_scaling import sha256
except ImportError:
    from rnn_scaling import sha256

def run(path):
    import numpy as np
    path=Path(path);cfg=json.loads(path.read_text());source=Path(cfg['input_csv'])
    if sha256(source)!=cfg['input_sha256']:raise ValueError('Wrong source table')
    with source.open() as f:rows=list(csv.DictReader(f))
    if len(rows)!=cfg['expected_rows']:raise ValueError('Unexpected row count')
    for r in rows:
        for k in ['accuracy','PQ_DC_mean','PS_DC_mean','severity']:r[k]=float(r[k])
        if not all(np.isfinite(r[k]) for k in ['accuracy','PQ_DC_mean','PS_DC_mean']):raise ValueError('Nonfinite row')
    out=Path(cfg['output_dir'])
    if out.exists() and any(out.iterdir()):raise FileExistsError('Use a new export directory')
    methods=cfg['methods']
    colors=[[0,.447,.741],[.85,.325,.098],[.466,.674,.188],[.494,.184,.556],[.929,.694,.125],[.301,.745,.933]]
    curves=[]
    for method,color in zip(methods,colors):
        levels=sorted({r['severity'] for r in rows if r['method']==method})
        for label,key,c in [('Accuracy','accuracy',color),('PQ','PQ_DC_mean',[.1,.1,.1]),('PS','PS_DC_mean',[.45,.45,.45])]:
            v=[np.array([r[key] for r in rows if r['method']==method and r['severity']==s]) for s in levels]
            if any(len(a)!=cfg['expected_rows_per_setting'] for a in v):raise ValueError('Unbalanced source settings')
            sem=[float(a.std(ddof=1)/np.sqrt(len(a))) for a in v]
            curves.append({'type':'export.ErrorBar','axes_title':method,'properties':dict(XData=list(range(1,len(levels)+1)),
                YData=[float(a.mean()) for a in v],YPositiveDelta=sem,YNegativeDelta=sem,DisplayName=label,Color=c),
                'severity_levels':levels,'sem_unit':'historical model-trial rows; not independent models'})
    baseline=[r for r in rows if [r['method'],r['severity']] in cfg['baseline_specs']]
    layers=[('baseline',baseline,[.466,.674,.188])]
    for layer in cfg['layers']:
        method,s,c=layer['method'],layer['severity'],layer['color']
        layers.append((method+' '+str(s),[r for r in rows if r['method']==method and r['severity']==s],c))
    scatter=[]
    for name,rr,color in layers:
        scatter.append({'type':'export.Scatter','properties':dict(XData=[r['PQ_DC_mean'] for r in rr],YData=[r['accuracy'] for r in rr],DisplayName=name,CData=color)})
    def quantile(key):return float(np.interp(cfg['baseline_quantile'],(np.arange(len(baseline))+.5)/len(baseline),np.sort([r[key] for r in baseline])))
    for axis,value in [('x',quantile('PQ_DC_mean')),('y',quantile('accuracy')-cfg['accuracy_tolerance'])]:
        scatter.append({'type':'export.ConstantLine','properties':{'Value':value,'Orientation':axis}})
    fits={'scope':'Descriptive Fig6g-h means, SEM and thresholds only; inference is supplied by the revised statistics module.'}
    out.mkdir(parents=True,exist_ok=True)
    for name,objects in [('method_parameter_trends_accuracy_PQ_PS',curves),('rnn_PQAcc_permutation_composer',scatter)]:
        (out/(name+'.json')).write_text(json.dumps({'source_status':cfg['source_status'],'objects':objects},indent=2))
    (out/'summary_scope.json').write_text(json.dumps(fits,indent=2))
    (out/'provenance.json').write_text(json.dumps({'source_sha256':sha256(source),'config_sha256':sha256(path),
        'code_sha256':sha256(__file__),'slurm_job_id':os.getenv('SLURM_JOB_ID'),'source_status':cfg['source_status']},indent=2))
    print(json.dumps(fits))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
