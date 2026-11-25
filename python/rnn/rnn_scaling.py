"""Export RNN scaling panel data on a cluster compute node; never render figures."""

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def reduce_scaling(source,sink):
    """Historical reduction: time/trial mean, D/C mean, pooled load z-score.

    Expected axes: time, trial, condition, model. Conditions are fixation,
    0-back and 2-back. Models remain paired across conditions. No hypothesis
    test is performed here; the historical unpaired test is not reproduced.
    """
    import numpy as np
    if source.shape!=sink.shape or source.ndim!=4 or source.shape[2]!=3:
        raise ValueError('Expected matching [time, trial, 3 conditions, model] arrays')
    if not np.isfinite(source).all() or not np.isfinite(sink).all():
        raise ValueError('Nonfinite inputs require an explicit exclusion policy')
    x=(source.mean(axis=(0,1))+sink.mean(axis=(0,1)))/2
    pooled=x[1:].reshape(-1)
    sd=pooled.std(ddof=1)
    if not sd>0:raise ValueError('Zero pooled variance')
    return (x[1:]-pooled.mean())/sd


def run(config_path):
    import numpy as np
    import scipy
    from scipy.io import loadmat
    config_path=Path(config_path).resolve()
    cfg=json.loads(config_path.read_text())
    root=Path(cfg['source_root']).resolve()
    out=Path(cfg['output_dir']).resolve()
    if out==root or root in out.parents:
        raise ValueError('Output must be outside the read-only source tree')
    if out.exists() and any(out.iterdir()):raise FileExistsError('Use a new empty export directory')
    # Verify the entire source contract before reading any arrays or writing output.
    checked=[]
    for item in cfg['inputs']:
        p=(root/item['path']).resolve()
        if root not in p.parents:raise ValueError('Source path escapes source_root')
        actual=sha256(p)
        if actual!=item['sha256']:raise ValueError('Source hash mismatch: '+item['path'])
        checked.append((item,p))
    exports=[]
    for item,p in checked:
        keys=('NumSource','NumSink') if item['metric']=='PQ' else ('SizeSource','SizeSink')
        d=loadmat(p,variable_names=keys)
        if list(d[keys[0]].shape)!=item['shape']:raise ValueError('Source shape mismatch')
        z=reduce_scaling(d[keys[0]],d[keys[1]])
        exports.append({'size':item['size'],'metric':item['metric'],'source_sha256':item['sha256'],
                        'n_models':z.shape[1],'z0back':z[0].tolist(),'z2back':z[1].tolist(),
                        'means':z.mean(axis=1).tolist(),'sem':(z.std(axis=1,ddof=1)/np.sqrt(z.shape[1])).tolist()})
    out.mkdir(parents=True,exist_ok=True)
    data=out/'scaling_plot_data.json'
    data.write_text(json.dumps({'schema_version':1,'status':cfg['source_status'],'panels':['F6e','F6f'],
        'boundary':'Descriptive historical normalization. No inferential statistics.', 'series':exports},indent=2))
    provenance={'config_sha256':sha256(config_path),'code_sha256':sha256(__file__),
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'slurm_job_id':os.environ.get('SLURM_JOB_ID', 'local'),'plot_data_sha256':sha256(data),
        'source_files':[{'path':i['path'],'sha256':i['sha256']} for i,_ in checked]}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2))
    print(json.dumps({'status':'exported','series':len(exports),'output_dir':str(out)}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
