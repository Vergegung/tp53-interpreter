"""打包已核验网站和六目录归档；不携带患者数据或开发缓存。"""
from pathlib import Path
import shutil
import hashlib
import json
import zipfile

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT.parent/'releases'/'TP53解读器_统一评估_v04_20261006'
if DEST.exists():
    raise FileExistsError('已有归档，避免覆盖，请使用新版本目录。')
DEST.mkdir(parents=True)
folders=['01_原始数据','02_清洗数据','03_图表输出','04_统计结果','05_分析报告','06_代码归档']
for d in folders:(DEST/d).mkdir()
for f in ['public.json.gz','manifest.json']:shutil.copy2(ROOT/'docs/data'/f,DEST/folders[0]/f)
for f in ['all_7467_hierarchical_v04.csv','representative_v04.csv','weight_sensitivity_v04.csv']:
    shutil.copy2(ROOT/'qa'/f,DEST/folders[1]/f)
for f in (ROOT/'qa/figures').iterdir():
    if f.is_file():shutil.copy2(f,DEST/folders[2]/f.name)
for f in ROOT.glob('qa/*.jpg'):shutil.copy2(f,DEST/folders[2]/f.name)
for f in ['index_audit_v04.json','node_tests_v04.txt','UI_QA_v04.md']:
    shutil.copy2(ROOT/'qa'/f,DEST/folders[3]/f)
shutil.copytree(ROOT/'qa/browser_exports',DEST/folders[3]/'browser_exports')
for source,name in [(ROOT/'REPORT_v04.md','优化与验证报告_v04.md'),(ROOT/'docs/METHODS.md','分层评估方法_v04.md'),(ROOT/'README.md','GitHub部署说明.md')]:
    shutil.copy2(source,DEST/folders[4]/name)
code=DEST/folders[5]/'github_site'
def excluded(directory,names):
    return set(names)&{'.git','.cache','__pycache__','public.json','.DS_Store','qa'}
shutil.copytree(ROOT,code,ignore=excluded)
# 仓库测试夹具必须完整；qa 是可重建的核验输出，归档在外层。
(code/'qa').mkdir()
shutil.copy2(ROOT/'qa/UI_QA_v04.md',code/'qa/UI_QA_v04.md')
shutil.copy2(ROOT/'qa/node_tests_v04.txt',code/'qa/node_tests_v04.txt')
readme='''# TP53 统一评估 v0.4 归档

浏览器版和核验已完成，真实 GitHub Pages 发布待登录。访问者无需本地执行；
后续发布使用 06_代码归档/github_site/，工作流文件在隐藏的 .github/workflows/ 下。

- 05_分析报告/优化与验证报告_v04.md：主报告、覆盖和科学边界。
- 03_图表输出/：实际页面预览和 600 DPI PNG、可编辑 PDF/SVG。
- 02_清洗数据/：7,467 个分层审计与权重敏感性。
- 04_统计结果/：软件和浏览器核验；本阶段没有临床统计检验。
- 患者 PDO/RNAseq/临床数据未打入归档；v0.3 与原始需求文档保留。

CHECKSUMS.sha256 覆盖所有归档文件（自身和 ARCHIVE_QC.json 除外）。
'''
(DEST/'开始阅读.md').write_text(readme)
zip_path=DEST/folders[5]/'tp53-interpreter-github-v04.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(code.rglob('*')):
        if p.is_file():z.write(p,p.relative_to(code))
hashes={}
for p in sorted(DEST.rglob('*')):
    if p.is_file():
        hashes[str(p.relative_to(DEST))]=hashlib.sha256(p.read_bytes()).hexdigest()
(DEST/'CHECKSUMS.sha256').write_text(''.join(v+'  '+k+'\n' for k,v in hashes.items()))
assert all(hashlib.sha256((DEST/k).read_bytes()).hexdigest()==v for k,v in hashes.items())
with zipfile.ZipFile(zip_path) as z:
    assert z.testzip() is None
    assert '.github/workflows/pages.yml' in z.namelist()
    assert 'docs/data/public.json.gz' in z.namelist()
    assert 'docs/data/public.json' not in z.namelist()
qc={'files_hash_verified':len(hashes),'zip_entries':len(z.namelist()),'patient_data_exported':False,
    'online_deployment_verified':False,'deployment_blocker':'GitHub browser sign-in and Pages setup pending',
    'snapshot_id':json.loads((ROOT/'qa/index_audit_v04.json').read_text())['snapshot_id']}
(DEST/'ARCHIVE_QC.json').write_text(json.dumps(qc,ensure_ascii=False,indent=2)+'\n')
print(DEST)
print(json.dumps(qc,ensure_ascii=False,indent=2))
