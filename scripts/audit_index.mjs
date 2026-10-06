// 冻结每个公开替换的分层值、覆盖与权重敏感性，便于回溯图表。
import fs from 'node:fs';
import zlib from 'node:zlib';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {interpret,DEFAULT_WEIGHTS} from '../docs/engine.mjs';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const db=JSON.parse(zlib.gunzipSync(fs.readFileSync(path.join(root,'docs/data/public.json.gz'))));
const out=path.join(root,'qa');fs.mkdirSync(out,{recursive:true});
const csv=(name,rows)=>{const keys=Object.keys(rows[0]);const esc=v=>v===null||v===undefined?'':JSON.stringify(String(v));fs.writeFileSync(path.join(out,name),keys.map(esc).join(',')+'\n'+rows.map(r=>keys.map(k=>esc(r[k])).join(',')).join('\n')+'\n');};
const summary={version:'0.4.0',snapshot_id:db.snapshot_id,weights:DEFAULT_WEIGHTS,variants:0,reportable:0,full_weighted_coverage:0,measured_structure_both_axes:0,status_counts:{},family_counts:{},coverage_counts:{}};
function row(q,w=DEFAULT_WEIGHTS){
 const r=interpret(q,db,{},w),h=r.hierarchical;
 return {query:q,protein:r.normalization.hgvs_p,kind:r.normalization.kind,structure_low:h.structure.coverage?h.structure.low:null,structure_high:h.structure.coverage?h.structure.high:null,structure_coverage:h.structure.coverage,function_low:h.function.coverage?h.function.low:null,function_high:h.function.coverage?h.function.high:null,function_coverage:h.function.coverage,dne_low:h.dne.coverage?h.dne.low:null,dne_high:h.dne.coverage?h.dne.high:null,dne_coverage:h.dne.coverage,total_compatible_low:h.total.low,total_compatible_high:h.total.high,coverage:h.total.coverage,functional_family_n:h.functional_family_n,reportable:h.reportable,status:h.status,applicability:h.allele_status,summary:r.unified.summary};
}
const all=Object.keys(db.variants).map(q=>row(q));
for(const r of all){summary.variants++;if(r.reportable)summary.reportable++;if(r.coverage>=1-1e-10)summary.full_weighted_coverage++;if(r.structure_coverage>=1-1e-10)summary.measured_structure_both_axes++;summary.status_counts[r.status]=(summary.status_counts[r.status]||0)+1;summary.family_counts[r.functional_family_n]=(summary.family_counts[r.functional_family_n]||0)+1;const c=r.coverage.toFixed(6);summary.coverage_counts[c]=(summary.coverage_counts[c]||0)+1;}
const examples=['R175H','R273H','R248Q','R280K','Y220C','G245S','R249S','R282W','C176S','C176F','P72R','c.586C>T','c.563del','c.376-2A>T'];
const representative=examples.map(q=>row(q));
const presets={function_priority:DEFAULT_WEIGHTS,equal:{structure:25,function:25,dne:25},structure_priority:{structure:50,function:25,dne:25}};
const sensitivity=examples.flatMap(q=>Object.entries(presets).map(([preset,w])=>({preset,w_structure:w.structure,w_function:w.function,w_dne:w.dne,...row(q,w)})));
csv('all_7467_hierarchical_v04.csv',all);csv('representative_v04.csv',representative);csv('weight_sensitivity_v04.csv',sensitivity);
fs.writeFileSync(path.join(out,'index_audit_v04.json'),JSON.stringify(summary,null,2)+'\n');
console.log(JSON.stringify(summary,null,2));
