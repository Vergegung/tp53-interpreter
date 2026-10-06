// 浏览器内计算：公开实验画像、条件性序列机制与透明的研究指数。
export const VERSION='0.4.0';
export const DEFAULT_WEIGHTS={structure:25,function:50,dne:25};
export const PROMOTERS=['CDKN1A','MDM2','BAX','SFN','TP53AIP1','GADD45A','NOXA','RRM2B'];
const aa3=['Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val','Ter'];
const aa1='ARNDCQEGHILKMFPSTWYV*';
const code={};let n=0;const acids='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
for(const a of 'TCAG')for(const b of 'TCAG')for(const c of 'TCAG')code[a+b+c]=acids[n++];
const translate=dna=>{let out='';for(let i=0;i+2<dna.length;i+=3)out+=code[dna.slice(i,i+3)];return out;};
const clone=x=>structuredClone(x);
const valid=x=>typeof x==='number'&&Number.isFinite(x);

export function normalize(input,db){
  let s=String(input).trim(),version='6',explicit=null;
  if(s.startsWith('NM_')){[explicit,s]=s.split(':');if(!['NM_000546.5','NM_000546.6'].includes(explicit)||!s)throw Error('仅支持 NM_000546.5/.6');version=explicit.slice(-1);}
  let p,kind,pos,extra={};const exact=s.startsWith('c.');
  if(exact){
    const snv=/^c\.(\d+)([ACGT])>([ACGT])$/.exec(s),indel=/^c\.(\d+)(?:_(\d+))?(del|dup|ins)([ACGT]*)$/.exec(s);
    if(snv){
      const base=+snv[1],ref=snv[2],alt=snv[3];
      if(base<1||base>1182||db.reference.cds[base-1]!==ref)throw Error('编码参考碱基与 RefSeq 不一致');
      pos=Math.floor((base-1)/3)+1;let codon=db.reference.cds.slice((pos-1)*3,pos*3).split('');codon[(base-1)%3]=alt;const mut=code[codon.join('')];
      if(pos===394){kind=mut==='*'?'synonymous':'stop_lost';p=mut==='*'?'p.*394=':`p.*394${mut}ext*?`;}
      else{const wild=db.reference.sequence[pos-1];kind=mut==='*'?'stop_gained':mut===wild?'synonymous':'missense';p=`p.${wild}${pos}${mut}`;if(pos===1&&mut!=='M')kind='start_lost';}
    }else if(indel){
      const a=+indel[1],b=+(indel[2]||indel[1]),op=indel[3],text=indel[4],tail=db.reference.tails[version];
      if(a<1||b<a||b>1179)throw Error('indel 当前支持 CDS 1–1179');
      if(op==='ins'&&(b!==a+1||!text))throw Error('ins 需相邻坐标与插入序列');
      if(op!=='ins'&&text&&tail.slice(a-1,b)!==text)throw Error('indel 参考序列不一致');
      const altered=op==='del'?tail.slice(0,a-1)+tail.slice(b):op==='dup'?tail.slice(0,b)+tail.slice(a-1,b)+tail.slice(b):tail.slice(0,a)+text+tail.slice(a);
      const translated=translate(altered),stop=translated.indexOf('*'),product=stop>=0?translated.slice(0,stop):translated;
      let j=0;for(;j<Math.min(product.length,393)&&product[j]===db.reference.sequence[j];j++);pos=j+1;
      if(pos>393){p=null;pos=null;kind='complex_coding';}
      else{
        const mut=product[pos-1]||'*',fs=(altered.length-tail.length)%3!==0;
        kind=mut==='*'?'stop_gained':fs?'frameshift':'inframe_indel';
        p=kind==='stop_gained'?`p.${db.reference.sequence[pos-1]}${pos}*`:kind==='frameshift'?`p.${db.reference.sequence[pos-1]}${pos}${mut}fs*${stop>=0?stop-pos+2:'?'}`:null;
        extra={ptc_c_position:stop>=0?stop*3+1:null,coding_length_delta:altered.length-tail.length,coding_edit_end:b,predicted_product_length:stop>=0?stop:null};
      }
    }else if(/^c\.(\d+)[+-]\d+[ACGT]>[ACGT]$/.test(s)){const base=+s.match(/^c\.(\d+)/)[1];if(base<1||base>1182)throw Error('剪接参考坐标超出支持范围');kind='splice_region';p=null;pos=null;}
    else throw Error('当前网页版未支持该复杂 HGVS；请使用明确的 SNV、del/dup/ins 或剪接区替换');
  }else{
    for(let i=0;i<aa3.length;i++)s=s.replaceAll(aa3[i],aa1[i]);s=s.replace(/^p\.\((.*)\)$/,'p.$1');if(!s.startsWith('p.'))s='p.'+s;
    const m=/^p\.([A-Z])(\d+)([A-Z*=])(fs(?:\*\d+|\*\?|\d+)?)?$/.exec(s);
    if(!m||+m[2]<1||+m[2]>393||db.reference.sequence[+m[2]-1]!==m[1]||!aa1.includes(m[3])&&m[3]!=='=')throw Error('蛋白 HGVS 无法解析或参考残基不一致');
    p=s;pos=+m[2];kind=m[4]?'frameshift':m[3]==='*'?'stop_gained':[m[1],'='].includes(m[3])?'synonymous':'missense';
  }
  const annotation=exact?db.coding_annotations[s]:null;
  let sourceConflict=false;
  if(exact&&/^c\.\d+[ACGT]>[ACGT]$/.test(s)&&pos<=393){
    for(const candidate of annotation?.protein||[]){try{const q=normalize(candidate,{...db,coding_annotations:{}});if(!(q.hgvs_p===p||q.kind===kind&&kind==='synonymous'))sourceConflict=true;}catch{sourceConflict=true;}}
  }
  return {input,hgvs_p:p||null,hgvs_c:exact?s:null,kind,position:pos,exact_nucleotide:exact,
    translation_transcript:exact?'NM_000546.'+version:null,transcript:exact?'NM_000546.'+version:null,
    input_transcript:explicit,domain:pos>=102&&pos<=292?'DBD':pos>=325&&pos<=356?'oligomerization':'other_or_unknown',
    source_protein_annotation_conflict:sourceConflict,...extra};
}

export function sequenceMechanism(norm,db,context={}){
  const k=norm.kind,pos=norm.position;let stop=norm.ptc_c_position??null;
  const r={variant_type:k,sequence_consequence:k,protein_presence:context.protein_presence||'unmeasured',
    ptc_c_position:null,ptc_exon:null,nmd:{status:'not_applicable',predicted:null},domains_if_translated:[],warnings:[],refolding_negative_control:'not_established'};
  if(!['unmeasured','present','absent'].includes(r.protein_presence))throw Error('蛋白状态需为未测、存在或缺失');
  if(r.protein_presence!=='unmeasured'){
    if(!context.protein_assay?.trim()||!context.protein_evidence_reference?.trim())throw Error('蛋白检测记录需要方法与来源');
    r.protein_evidence={assay:context.protein_assay,reference:context.protein_evidence_reference,provenance:'user_supplied_unverified'};
    if(r.protein_presence==='absent')r.refolding_negative_control='candidate_after_independent_protein_validation';
  }
  if(['frameshift','stop_gained'].includes(k)){
    if(k==='stop_gained'&&stop===null&&norm.exact_nucleotide)stop=(pos-1)*3+1;
    r.ptc_c_position=stop;
    if(stop===null)r.nmd={status:'unresolved_ptc_or_transcript',predicted:null};
    else{
      const version=norm.translation_transcript?.slice(-1)||'6',exons=db.reference.exons[version],last=exons.at(-2).end;
      const delta=norm.coding_length_delta||0,junction=(norm.coding_edit_end||0)<=last?last+delta:last,distance=junction-(stop+2);
      const predicted=distance>55?true:distance<50?false:null;
      r.ptc_exon=exons.find(e=>e.start<=stop-delta&&stop-delta<=e.end)?.number??null;
      r.nmd={status:stop<=450?'positional_rule_with_reinitiation_caveat':predicted===null?'boundary_50_55_nt':'positional_prediction_only',predicted,stop_end_to_last_junction_nt:distance,last_junction_c_position:junction,transcript:'NM_000546.'+version};
      if(stop<=450)r.warnings.push('早期 PTC 可能有翻译再起始或亚型表达；位置预测不等于蛋白缺失。');
    }
    for(const [domain,a,b] of [['DBD',102,292],['oligomerization',325,356]]){
      const end=stop!==null?Math.min(Math.floor((stop-1)/3),pos-1):(pos?pos-1:0);
      r.domains_if_translated.push({domain,status:k==='frameshift'&&pos<=b?'lost_or_altered':end<a?'lost':end<b?'partial':'sequence_retained'});
    }
    r.warnings.push('截断的序列后果可计算；实际 RNA、蛋白及亚型需验证，不能按外显子位置判定零蛋白。');
  }else if(k==='splice_region'){r.sequence_consequence='splice_region_variant_RNA_effect_unresolved';r.nmd={status:'requires_observed_RNA_product',predicted:null};r.warnings.push('实际 RNA 产物未决，不预设移码或 NMD。');}
  else if(k==='inframe_indel'){r.sequence_consequence='reading_frame_preserved_local_sequence_changed';r.warnings.push('框内局部序列改变；不能借用错义实验分值。');}
  else if(k==='start_lost')r.warnings.push('起始丢失需确认替代起始与亚型，不使用 M1 普通错义分值。');
  else if(k==='stop_lost')r.warnings.push('终止丢失可能产生延长蛋白，停止位置未定，不作为截断。');
  else if(k==='synonymous')r.warnings.push('蛋白序列不变，仍需排除剪接与表达作用。');
  return r;
}

export function axisInterval(rank){
  if(valid(rank))return {low:rank,high:rank,coverage:1,point:rank};
  if(Array.isArray(rank)&&rank.length===2&&rank.every(valid))return {low:Math.min(...rank),high:Math.max(...rank),coverage:1,point:rank[0]===rank[1]?rank[0]:null};
  return {low:0,high:100,coverage:0,point:null};
}
export function aggregate(components){
  const total=components.reduce((a,c)=>a+c.weight,0);if(total<=0)throw Error('权重总和必须大于零');
  let low=0,high=0,covered=0,observedLow=0,observedHigh=0;
  for(const c of components){const w=c.weight/total,v=c.interval;low+=w*v.low;high+=w*v.high;covered+=w*v.coverage;
    observedLow+=w*(v.low-(1-v.coverage)*0);observedHigh+=w*(v.high-(1-v.coverage)*100);}
  return {low,high,coverage:covered,point:covered===1&&Math.abs(high-low)<1e-9?low:null,
    observed_low:covered>0?observedLow/covered:null,observed_high:covered>0?observedHigh/covered:null,
    interval_kind:'fixed-weight compatibility range; not CI',components};
}
export function scoreProfile(ranks,weights=DEFAULT_WEIGHTS){
  if(!Object.keys(DEFAULT_WEIGHTS).every(k=>valid(weights[k])&&weights[k]>=0)||Object.values(weights).reduce((a,b)=>a+b,0)<=0)throw Error('分项权重必须非负且总和大于零');
  const axis=ranks.map(axisInterval),family=(name,indices)=>({name,interval:aggregate(indices.map(i=>({name:String(i),weight:1,interval:axis[i]})))});
  const structure=aggregate([0,1].map(i=>({name:i===0?'apo 稳定性':'锌亲和力',weight:1,interval:axis[i]})));
  const families=[family('Kato',[2]),family('Giacomelli LOF',[3,4]),family('Kotler',[6])];
  const functional=aggregate(families.map(f=>({...f,weight:1})));const dne=aggregate([{name:'Giacomelli WT/nutlin',weight:1,interval:axis[5]}]);
  const total=aggregate([{name:'结构',weight:weights.structure,interval:structure},{name:'功能缺失',weight:weights.function,interval:functional},{name:'显性负性',weight:weights.dne,interval:dne}]);
  const functionalFamilyN=families.filter(f=>f.interval.coverage>0).length;
  return {formula_version:'research-index-0.4.0',weights:{...weights},structure,function: functional,dne,total,
    functional_family_n:functionalFamilyN,reportable:total.coverage>=.6&&functionalFamilyN>=2,
    status:total.coverage===0?'not_quantifiable':total.coverage<.6||functionalFamilyN<2?'insufficient_coverage':'descriptive_research_index',
    note:'预设权重的描述性研究指数，尚未外部校准。区间包含缺失与来源歧义，不是概率或置信区间。'};
}

export function interpret(input,db,context={},weights=DEFAULT_WEIGHTS){
  const norm=normalize(input,db),record=norm.kind==='missense'?db.variants[norm.hgvs_p?.slice(2)]:null;
  const classes=record?Object.fromEntries(db.class_fields.map((k,i)=>[k,record.classes[i]])):{};
  const ranks=record?.ranks||Array(7).fill(null),non=sequenceMechanism(norm,db,context);
  const annotation=norm.exact_nucleotide?db.coding_annotations[norm.hgvs_c]:null;
  const spliceMax=context.spliceai_max??annotation?.spliceai_max??null;
  if(spliceMax!==null&&(!valid(spliceMax)||spliceMax<0||spliceMax>1))throw Error('SpliceAI 需在 0–1 之间');
  const blocked=norm.source_protein_annotation_conflict||valid(spliceMax)&&spliceMax>=.2||context.rna_splicing_abnormal===true;
  const unresolved=norm.exact_nucleotide&&context.splicing_ruled_out!==true&&!annotation?.spliceai_complete;
  const score=scoreProfile(ranks,weights);
  if(blocked){score.allele_status='hold_identity_or_splicing';score.allele_score=null;}
  else score.allele_status=norm.exact_nucleotide?(unresolved?'coding_identity_splicing_unresolved':'coding_identity'): 'protein_substitution_only';
  if(blocked&&['stop_gained','frameshift'].includes(norm.kind)){non.nmd.canonical_splicing_assumption_unverified=true;non.warnings.push('存在身份或剪接门控；PTC/NMD 为 canonical 剪接保持时的假设后果。');}
  const manual=context.dimension_inputs||[];if(!Array.isArray(manual)||manual.length>20)throw Error('补充测量最多 20 条');
  const checked=manual.map(m=>{if(!['structure','function'].includes(m.dimension)||!['metric','unit','condition','source'].every(k=>typeof m[k]==='string'&&m[k].trim())||!valid(m.value))throw Error('补充测量需完整的维度、指标、有限数值、单位、条件、来源');return {...m,provenance:'user_supplied_unverified'};});
  const kindText={stop_gained:'提前终止：RNA 降解或截短蛋白候选',frameshift:'阅读框改变：PTC 与蛋白产物需结合验证',splice_region:'剪接机制待 RNA 产物确认',start_lost:'canonical 起始丢失：替代起始待确认',stop_lost:'正常终止丢失：C 端延长候选',inframe_indel:'框内局部序列改变：功能未定',synonymous:'canonical 序列保留：剪接与表达作用待排除',complex_coding:'canonical 蛋白末端之后的后果未解析'};
  const labels={LOF:'Kato 与 Giacomelli 支持功能缺失',retained:'Kato 与 Giacomelli 支持功能保留',insufficient_assays:'直接功能实验覆盖不足',discordant_or_partial:'部分功能或来源分歧',unassessed:'缺乏直接功能实验'};
  const summary=blocked?'核苷酸身份或剪接存在冲突；蛋白画像仅作参考':kindText[norm.kind]||labels[classes.consensus_class]||'当前功能未定';
  const evidence=blocked?'需解决身份 / RNA':norm.kind!=='missense'?'序列后果 + 条件性机制':score.functional_family_n>=2?'多实验家族支持':score.functional_family_n===1?'单实验家族':'计算或未覆盖';
  return {version:VERSION,snapshot_id:db.snapshot_id,normalization:norm,classes,promoters:record?.promoters||Array(8).fill(null),
    raw:record?.raw||Array(7).fill(null),ranks,structure:record?db.structures[norm.hgvs_p.slice(2)]||null:null,
    kotler_records:record?.kotler_records||[],conflicts:record?.conflicts||[],non_missense:non,
    hierarchical:score,unified:{summary,evidence,mechanism_route:norm.kind==='missense'?'experiment_profile':'sequence_mechanism',
      applicability:score.allele_status,clinical_response_probability:null},manual:checked};
}
