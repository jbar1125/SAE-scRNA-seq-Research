#!/usr/bin/env python3
"""Typeset the existing internal draft; preserve all author/disclosure placeholders."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out-dir', required=True, type=Path)
    args = p.parse_args()
    if args.out_dir.exists():
        p.error('choose a fresh output directory')
    for tool in ['pandoc', 'xelatex', 'pdftotext', 'pdfinfo']:
        if shutil.which(tool) is None:
            p.error(f'{tool} is required for this optional review-artifact workflow')
    out = args.out_dir.resolve()
    build = out / 'build'
    build.mkdir(parents=True)
    source = ROOT / 'docs/MANUSCRIPT_DRAFT.md'
    fingerprint = sha(source)
    lines = source.read_text().splitlines()
    figures = [(ROOT / 'docs' / name).resolve() for name in
               re.findall(r'!\[[^\]]*\]\(([^)]+)\)', '\n'.join(lines))]
    figure_hashes = {str(path.relative_to(ROOT)): sha(path) for path in figures}
    title = lines[0].removeprefix('# ')
    body = build / 'manuscript-body.md'
    body.write_text('\n'.join(lines[1:]) + '\n')
    layout = build / 'review-layout.lua'
    layout.write_text('''function Header(el)
  local title = pandoc.utils.stringify(el.content)
  if title == "Abstract" or title == "Supplementary figure" then
    return {pandoc.RawBlock("latex", "\\\\clearpage"),
      pandoc.RawBlock("openxml", '<w:p><w:r><w:br w:type="page"/></w:r></w:p>'), el}
  end
end
function Image(el)
  if el.src:match("fig_robustness%.png$") then el.attributes.width = "6in" end
  return el
end
''')
    env = os.environ.copy()
    for name, directory in [('XDG_CACHE_HOME', 'cache'), ('TEXMFVAR', 'texmf-var'), ('TEXMFCONFIG', 'texmf-config')]:
        path = build / directory
        path.mkdir()
        env[name] = str(path)
    base = ['pandoc', str(body), '--from=markdown-implicit_figures', '--standalone',
            '--lua-filter=' + str(layout),
            '--shift-heading-level-by=-1', '--resource-path=' + str(ROOT / 'docs'),
            '--metadata', 'title=' + title, '--metadata',
            'date=Internal review draft, ' + datetime.now(timezone.utc).date().isoformat()]
    with (build / 'render.log').open('w') as log:
        subprocess.run(base + ['--output', str(out / 'manuscript-review.docx')],
                       stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
        subprocess.run(base + ['--pdf-engine=xelatex', '--variable', 'geometry:margin=1in',
            '--variable', 'fontsize=11pt', '--variable', 'linestretch=1.15',
            '--variable', 'mainfont=DejaVu Serif', '--variable', 'sansfont=DejaVu Sans',
            '--variable', 'monofont=DejaVu Sans Mono', '--output', str(out / 'manuscript-review.pdf')],
            stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    if sha(source) != fingerprint:
        raise RuntimeError('manuscript changed during rendering; preserve outputs but rebuild from one version')
    if figure_hashes != {str(path.relative_to(ROOT)): sha(path) for path in figures}:
        raise RuntimeError('a figure changed during rendering; preserve outputs and rebuild')
    text = subprocess.check_output(['pdftotext', str(out / 'manuscript-review.pdf'), '-'], text=True)
    compact = ''.join(text.split())
    for fragment in ['Abstract', '428,505', '105', 'Contributions,funding,anddisclosures', 'References']:
        if fragment not in compact:
            raise RuntimeError(f'rendered PDF is missing expected content: {fragment}')
    with zipfile.ZipFile(out / 'manuscript-review.docx') as package:
        if package.testzip() is not None:
            raise RuntimeError('the generated Word package failed its CRC check')
        image_count = sum(n.startswith('word/media/') for n in package.namelist())
        if image_count != 5:
            raise RuntimeError('the five manuscript figures were not embedded')
    outputs = {p.name: {'bytes': p.stat().st_size, 'sha256': sha(p)}
               for p in [out / 'manuscript-review.pdf', out / 'manuscript-review.docx']}
    manifest = {'source': str(source.relative_to(ROOT)), 'source_sha256': fingerprint,
        'generator_sha256': sha(Path(__file__).resolve()),
        'layout_filter_sha256': sha(layout),
        'figure_sha256': figure_hashes,
        'rendered_utc': datetime.now(timezone.utc).isoformat(), 'outputs': outputs,
        'embedded_figures': image_count, 'pdf_content_checks': 'PASS',
        'pandoc_version': subprocess.check_output(['pandoc', '--version'], text=True).splitlines()[0],
        'pdf_info': subprocess.check_output(['pdfinfo', str(out / 'manuscript-review.pdf')], text=True),
        'scope': 'Internal review copies. Placeholders remain; no author approval, submission or venue formatting is implied.'}
    (out / 'render_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'outputs': outputs, 'embedded_figures': image_count}), flush=True)


if __name__ == '__main__':
    main()
