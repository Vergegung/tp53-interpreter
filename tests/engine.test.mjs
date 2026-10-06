import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import crypto from 'node:crypto';
import zlib from 'node:zlib';
import {normalize,interpret,axisInterval,aggregate,scoreProfile} from '../docs/engine.mjs';
const file=new URL('../docs/data/public.json.gz',import.meta.url),raw=zlib.gunzipSync(fs.readFileSync(file)),data=JSON.parse(raw.toString('utf8'));
function csv(file){const text=fs.readFileSync(file,'utf8').trim(),rows=[];let row=[],v='',quoted=false;for(let i=0;i<text.length;i++){const c=text[i];if(c==='"'){if(quoted&&text[i+1]==='"'){v+='"';i++;}else quoted=!quoted;}else if(c===','&&!quoted){row.push(v);v='';}else if(c==='\n'&&!quoted){row.push(v.replace(/\r$/,''));rows.push(row);row=[];v='';}else v+=c;}row.push(v.replace(/\r$/,''));rows.push(row);const headers=rows.shift();return rows.map(r=>Object.fromEntries(headers.map((h,i)=>[h,r[i]])));}
const approx=(a,b)=>assert.ok(Math.abs(a-b)<1e-8,`${a} != ${b}`);

test('public data matches SHA256 and has exact variant reference',()=>{
 const manifest=JSON.parse(fs.readFileSync(new URL('../docs/data/manifest.json',import.meta.url)));
 assert.equal(crypto.createHash('sha256').update(raw).digest('hex'),manifest.public_json_sha256);
 assert.equal(Object.keys(data.variants).length,7467);assert.equal(manifest.structural_measured_n,23);assert.equal(manifest.patient_data_exported,false);
 assert.deepEqual(data.axis_specs.map(v=>v.reference_n),[23,20,2314,7448,7448,7448,4069]);
 assert.ok(data.axis_specs.every(v=>v.unit));
 for(const p of Object.keys(data.variants)){const n=normalize(p,data);assert.equal(n.kind,'missense');for(const a of data.variants[p].ranks){if(a!==null)for(const v of Array.isArray(a)?a:[a])assert.ok(v>=0&&v<=100);}}
});
test('all 3546 coding SNVs agree with the audited Python core',()=>{
 const rows=csv(new URL('fixtures/all_coding_snv_v03.csv',import.meta.url));assert.equal(rows.length,3546);
 for(const r of rows){const n=normalize(r.hgvs_c,data);assert.equal(n.hgvs_p,r.observed_hgvs_p,r.hgvs_c);assert.equal(n.kind,r.variant_type,r.hgvs_c);assert.equal(n.source_protein_annotation_conflict,r.historical_conflict==='True',r.hgvs_c);}
});
test('both RefSeq tails: all 4716 del/dup consequences and PTCs agree',()=>{
 const rows=csv(new URL('fixtures/all_single_base_indels_v03.csv',import.meta.url));assert.equal(rows.length,4716);
 for(const r of rows){const n=normalize(r.hgvs_c,data);assert.equal(n.hgvs_p,r.hgvs_p||null,r.hgvs_c);assert.equal(n.kind,r.kind,r.hgvs_c);assert.equal(n.translation_transcript,r.translation_transcript);if(r.interpreter_stop_c_position)assert.equal(n.ptc_c_position,Number(r.interpreter_stop_c_position),r.hgvs_c);}
});
test('representative source axes and NMD retain canonical core behavior',()=>{
 const fixtures=JSON.parse(fs.readFileSync(new URL('fixtures/representative.json',import.meta.url)));
 for(const v of fixtures){const query=v.normalization.hgvs_c?`${v.normalization.transcript}:${v.normalization.hgvs_c}`:v.normalization.hgvs_p;
  if(!query)continue;const r=interpret(query,data);assert.equal(r.normalization.kind,v.normalization.kind,query);
  assert.equal(r.non_missense.nmd.predicted,v.non_missense.nmd.predicted,query);
  assert.equal(r.non_missense.ptc_c_position,v.non_missense.ptc_c_position,query);
  for(let i=0;i<7;i++){const a=v.profile.axes[i],expected=a.percentile??a.percentile_range;if(Array.isArray(expected)){assert.deepEqual(r.ranks[i],expected);}else assert.equal(r.ranks[i],expected??null);}
 }
});
test('hierarchical weights do not count 8 promoters or 2 Gia conditions as extra families',()=>{
 const s=scoreProfile([20,40,60,70,90,50,[30,50]]);approx(s.structure.low,30);approx(s.function.low,170/3);approx(s.function.high,190/3);
 approx(s.total.low,.25*30+.5*170/3+.25*50);approx(s.total.high,.25*30+.5*190/3+.25*50);
 assert.equal(s.functional_family_n,3);assert.equal(s.total.point,null);
});
test('missing observations expand compatibility range without zero imputation',()=>{
 const s=scoreProfile([null,null,80,80,80,80,80]);approx(s.total.low,60);approx(s.total.high,85);approx(s.total.coverage,.75);
 approx(s.total.observed_low,80);approx(s.total.observed_high,80);assert.ok(s.reportable);
 const empty=scoreProfile(Array(7).fill(null));assert.equal(empty.status,'not_quantifiable');assert.equal(empty.total.coverage,0);assert.equal(empty.total.low,0);assert.equal(empty.total.high,100);assert.equal(empty.reportable,false);
});
test('all possible outcomes stay in 0-100 with valid nonnegative weights',()=>{
 for(let i=0;i<100;i++){const ranks=Array.from({length:7},(_,j)=>(i+j)%4===0?null:(i*37+j*11)%101);const s=scoreProfile(ranks,{structure:i%50,function:50,dne:25});assert.ok(s.total.low>=0&&s.total.high<=100+1e-10);assert.ok(s.total.low<=s.total.high);}
 assert.throws(()=>scoreProfile(Array(7).fill(1),{structure:0,function:0,dne:0}));assert.throws(()=>scoreProfile(Array(7).fill(1),{structure:-1,function:50,dne:25}));
});
test('changing valid weights changes research index and preserves raw data',()=>{
 const a=interpret('R175H',data,{}, {structure:100,function:0,dne:0}),b=interpret('R175H',data,{}, {structure:0,function:100,dne:0});
 assert.notEqual(a.hierarchical.total.low,b.hierarchical.total.low);assert.deepEqual(a.raw,b.raw);
});
test('non-missense does not borrow measurements or automatically imply no protein',()=>{
 for(const q of ['c.1A>G','c.1180T>C','c.563del','c.586C>T','c.4_6del','c.376-2A>T']){const r=interpret(q,data);assert.ok(r.ranks.every(x=>x===null));assert.equal(r.hierarchical.reportable,false);assert.equal(r.non_missense.protein_presence,'unmeasured');assert.equal(r.unified.clinical_response_probability,null);}
 assert.equal(interpret('c.1123C>T',data).non_missense.nmd.predicted,false);assert.equal(interpret('R196*',data).non_missense.nmd.predicted,null);
});
test('variant-specific structure not borrowed and censored DNA remains a bound',()=>{
 assert.equal(interpret('R282W',data).structure.measurements.length,0);
 const y=interpret('Y220C',data);assert.equal(y.ranks[1],null);approx(y.structure.measurements.find(m=>m.metric==='apo_destabilization_vs_wt').value,2.1);
 assert.ok(interpret('R273H',data).structure.dna_binding.every(m=>m.value===null&&m.lower_bound===25));
});
test('identity and splice gates hold allele-level application',()=>{
 assert.equal(interpret('c.216C>A',data).hierarchical.allele_status,'hold_identity_or_splicing');
 assert.equal(interpret('c.524G>A',data,{spliceai_max:.9}).hierarchical.allele_status,'hold_identity_or_splicing');
 assert.equal(interpret('R175H',data).hierarchical.allele_status,'protein_substitution_only');
});
test('user protein evidence and experimental values do not overwrite public scores',()=>{
 assert.throws(()=>interpret('R196*',data,{protein_presence:'absent'}));
 const r=interpret('Y220C',data,{dimension_inputs:[{dimension:'structure',metric:'delta_Tm',value:-5,unit:'C',condition:'DSF',source:'lab-001'}]});
 approx(r.raw[0],2.1);assert.equal(r.manual[0].provenance,'user_supplied_unverified');
 assert.throws(()=>interpret('Y220C',data,{dimension_inputs:[{value:4}]}));
});
test('invalid and unsupported HGVS do not silently switch identities',()=>{
 for(const q of ['R175Qbad','R176H','NM_000546.4:c.524G>A','c.524A>G','c.1_2ins','c.1_2000del'])assert.throws(()=>normalize(q,data),q);
});
