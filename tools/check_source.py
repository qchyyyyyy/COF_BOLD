"""Validate the public source inventory without importing scientific modules."""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_SUFFIXES = {'.mat', '.npy', '.npz', '.h5', '.hdf5', '.nii', '.gii',
                      '.fif', '.pt', '.pth', '.xlsx', '.docx', '.pdf', '.png',
                      '.jpg', '.tiff', '.csv', '.tsv', '.mexw64', '.mexa64', '.exe', '.dll'}
SKIP_PARTS = {'.git', '__pycache__', '.pytest_cache', '.venv'}


def files(root=ROOT):
    return [p for p in sorted(root.rglob('*')) if p.is_file()
            and not any(s in SKIP_PARTS for s in p.relative_to(root).parts)]


def check(root=ROOT, verify_manifest=False):
    errors = []
    paths = files(root)
    for p in paths:
        rel = p.relative_to(root).as_posix()
        if any(part in {'pipelines','experiments'} for part in p.relative_to(root).parts):
            errors.append(f'internal project layout leaked into release: {rel}')
        if not rel.startswith('third_party/') and re.search(r'P\d{2}-EXP-\d+', rel):
            errors.append(f'internal experiment ID in public path: {rel}')
        if p.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f'data or binary payload: {rel}')
        if '.local.' in p.name or p.name in {'.env', 'id_rsa', 'id_ed25519'}:
            errors.append(f'private configuration: {rel}')
        try:
            s = p.read_text(encoding='utf-8-sig')
        except UnicodeDecodeError:
            errors.append(f'non-text source file: {rel}')
            continue
        if p.suffix == '.py':
            try:
                ast.parse(s, filename=rel)
            except SyntaxError as exc:
                errors.append(f'Python syntax: {rel}: {exc}')
        # These rules identify embedded historical installation paths and
        # common credential forms; they are not a universal privacy proof.
        for pattern in (r'/public/(?:home|mig)', r'(?i)[CDE]:[\\/](?:Users|spatiotemporal_patterns|qcy)',
                        r'gh[pousr]_[A-Za-z0-9]{30,}', r'AKIA[0-9A-Z]{16}',
                        r'-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----'):
            if re.search(pattern, s):
                errors.append(f'private path or credential pattern: {rel}')
    names = {}
    for p in (root / 'matlab').rglob('*.m'):
        if p.stem in names:
            errors.append(f'duplicate MATLAB public function filename: {p.stem}')
        names[p.stem] = p
        if not re.fullmatch(r'[a-z][a-z0-9_]*', p.stem):
            errors.append(f'inconsistent MATLAB filename: {p.name}')
    if verify_manifest:
        m = root / 'SOURCE_MANIFEST.json'
        if not m.exists():
            errors.append('SOURCE_MANIFEST.json missing')
        else:
            manifest = json.loads(m.read_text(encoding='utf-8'))
            expected = {x['path']: x for x in manifest['files']}
            actual = {p.relative_to(root).as_posix(): p for p in paths
                      if p.name != 'SOURCE_MANIFEST.json'}
            if set(expected) != set(actual):
                errors.append(f'manifest file set mismatch: {sorted(set(expected)^set(actual))}')
            for name in expected.keys() & actual.keys():
                if hashlib.sha256(actual[name].read_bytes()).hexdigest() != expected[name]['sha256']:
                    errors.append(f'hash mismatch: {name}')
    return {'status': 'passed' if not errors else 'failed', 'file_count': len(paths),
            'python_files': sum(p.suffix == '.py' for p in paths), 'errors': errors}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', action='store_true')
    a = p.parse_args()
    report = check(verify_manifest=a.manifest)
    print(json.dumps(report, indent=2))
    sys.exit(bool(report['errors']))
