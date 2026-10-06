import {interpret,scoreProfile,DEFAULT_WEIGHTS,PROMOTERS} from './engine.mjs';
const $=id=>document.getElementById(id),colors=['#E64B35','#4DBBD5','#00A087','#3C5488','#7E6163','#F39B7F'];
let db,results=[],queries=[],current=0,weights={...DEFAULT_WEIGHTS},lastContext={};
const num=x=>x===null||x===undefined?'—':typeof x==='number'?(Math.abs(x)>0&&Math.abs(x)<.001?x.toExponential(3):Number(x.toPrecision(5)).toString()):String(x);
const pct=x=>Math.round(x*100)+'%';
const bounds=v=>Math.abs(v.high-v.low)<.05?v.low.toFixed(1):v.low.toFixed(1)+'–'+v.high.toFixed(1);
const rank=v=>Array.isArray(v)?v.map(num).join('–'):num(v);
const terms={missense:'错义',stop_gained:'无义 / 提前终止',frameshift:'移码',splice_region:'剪接区',inframe_indel:'框内插入 / 缺失',start_lost:'起始丢失',stop_lost:'终止丢失',synonymous:'同义',complex_coding:'复杂编码变异',protein_substitution_only:'蛋白替换画像；未指定核苷酸等位基因',coding_identity:'编码身份已核对；剪接影响仍需结合实验',coding_identity_splicing_unresolved:'编码身份已核对；剪接证据不完整',hold_identity_or_splicing:'身份或剪接冲突待解决；仅显示假设蛋白画像',unmeasured:'未测量',present:'检测到蛋白（用户提供，未核验）',absent:'未检测到蛋白（用户提供，未核验）',not_applicable:'不适用',not_established:'尚未建立',unresolved_ptc_or_transcript:'PTC 或转录本未明确',positional_rule_with_reinitiation_caveat:'位置预测；需考虑翻译再起始',boundary_50_55_nt:'距末次连接点 50–55 nt，位置规则未决',positional_prediction_only:'仅依据 PTC 位置预测',requires_observed_RNA_product:'需明确实际 RNA 产物',candidate_after_independent_protein_validation:'独立验证蛋白缺失后，可考虑作为候选',splice_region_variant_RNA_effect_unresolved:'剪接区变异，实际 RNA 后果未决',reading_frame_preserved_local_sequence_changed:'阅读框保留，局部序列改变',lost_or_altered:'原始序列丢失或改变',lost:'原始序列丢失',partial:'部分原始序列保留',sequence_retained:'原始序列保留',oligomerization:'寡聚化结构域'};
const readable=x=>x==='structure'?'结构':x==='function'?'功能':terms[x]||num(x);
const levelValue=v=>v.coverage>0?bounds(v):'NA（缺乏可比测量）';
function el(tag,text,cls,parent){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;if(parent)parent.append(n);return n;}
function table(parent,heads,rows){const wrap=el('div',undefined,'tablewrap',parent),t=el('table',undefined,undefined,wrap),h=el('tr',undefined,undefined,el('thead',undefined,undefined,t));heads.forEach(v=>el('th',v,undefined,h));const body=el('tbody',undefined,undefined,t);for(const row of rows){const tr=el('tr',undefined,undefined,body);row.forEach(v=>el('td',Array.isArray(v)?rank(v):num(v),undefined,tr));}}
function ref(parent,label,url){if(!/^https?:\/\//.test(url||''))return;const a=el('a',label,'ref',parent);a.href=url;a.target='_blank';a.rel='noopener noreferrer';}
function notice(parent,text){el('div',text,'notice',parent);}
function save(text,type,name){const u=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=u;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(u),10000);}
function svg(tag,attrs,text,parent){const n=document.createElementNS('http://www.w3.org/2000/svg',tag);Object.entries(attrs).forEach(([k,v])=>n.setAttribute(k,v));if(text!==undefined)n.textContent=text;parent?.append(n);return n;}
function chart(target,axes,series,max=100,wt=false){
 target.replaceChildren();const s=svg('svg',{viewBox:'0 0 660 490',xmlns:'http://www.w3.org/2000/svg',role:'img','aria-label':wt?'Kato eight promoter activity':'TP53 seven-axis source profile'},undefined,target);
 const cx=330,cy=230,r=144,n=axes.length,point=(i,v)=>[cx+Math.sin(i*2*Math.PI/n)*r*v/max,cy-Math.cos(i*2*Math.PI/n)*r*v/max],ps=a=>a.map(x=>x.join(',')).join(' ');
 for(let tick=1;tick<=4;tick++){const v=max*tick/4;svg('polygon',{points:ps(axes.map((_,i)=>point(i,v))),fill:'none',stroke:'#dce5e6','stroke-width':1},undefined,s);svg('text',{x:336,y:230-r*v/max+12,fill:'#86969e','font-family':'Arial','font-size':11},num(v),s);}
 if(wt)svg('polygon',{points:ps(axes.map((_,i)=>point(i,100))),fill:'none',stroke:'#899ca3','stroke-dasharray':'5 4','stroke-width':1.4},undefined,s);
 axes.forEach((a,i)=>{const p=point(i,max),t=point(i,max*1.28),anchor=Math.abs(t[0]-cx)<8?'middle':t[0]>cx?'start':'end';svg('line',{x1:cx,y1:cy,x2:p[0],y2:p[1],stroke:'#e2e9eb'},undefined,s);const lines=a.label.split(' / ');lines.forEach((l,j)=>svg('text',{x:t[0],y:t[1]+j*14,'text-anchor':anchor,'font-family':'Arial','font-size':12,fill:'#294550'},l,s));if(a.reference_n)svg('text',{x:t[0],y:t[1]+lines.length*14,'text-anchor':anchor,'font-family':'Arial','font-size':10,fill:'#84949a'},'N='+a.reference_n,s);});
 series.forEach((v,k)=>v.values.forEach((x,i)=>{const next=v.values[(i+1)%n];if(typeof x==='number'&&typeof next==='number'){const a=point(i,x),b=point((i+1)%n,next);svg('line',{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:colors[k],'stroke-width':1.8},undefined,s);}if(typeof x==='number'){const a=point(i,x),dot=svg('circle',{cx:a[0],cy:a[1],r:3.5,fill:colors[k],stroke:'white','stroke-width':.7},undefined,s);svg('title',{},v.name+' · '+axes[i].label+': '+num(x),dot);}if(Array.isArray(x)){const a=point(i,x[0]),b=point(i,x[1]);svg('line',{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:colors[k],'stroke-width':5,'stroke-opacity':.5},undefined,s);}}));
 axes.forEach((_,i)=>{if(series.every(v=>v.values[i]===null)){const a=point(i,max*.5);svg('text',{x:a[0],y:a[1],'text-anchor':'middle','font-family':'Arial','font-size':11,fill:'#91a0a7'},'NA',s);}});
 svg('text',{x:330,y:466,'text-anchor':'middle','font-family':'Arial','font-size':11,fill:'#6e858e'},wt?'Dashed line: WT = 100%; raw source activity':'Gaps = unavailable; bars = source record ranges',s);
}
function showCharts(){
 chart($('radar'),db.axis_specs,results.map(r=>({name:r.normalization.hgvs_p||r.normalization.hgvs_c,values:r.ranks})));
 const max=Math.max(150,Math.ceil(Math.max(100,...results.flatMap(r=>r.promoters).filter(x=>x!==null))/50)*50);
 chart($('promoters'),PROMOTERS.map(label=>({label})),results.map(r=>({name:r.normalization.hgvs_p||r.normalization.hgvs_c,values:r.promoters})),max,true);
 $('legend').replaceChildren();results.forEach((r,i)=>{const s=el('span',undefined,undefined,$('legend'));const dot=el('i',undefined,undefined,s);dot.style.background=colors[i];s.append(document.createTextNode(r.normalization.hgvs_p||r.normalization.hgvs_c));});
}
function levelCard(parent,title,v){const c=el('article',undefined,'panel levelcard',parent);el('h3',title,undefined,c);el('strong',v.coverage>0?bounds(v):'NA',undefined,c);el('p','加权覆盖 '+pct(v.coverage),undefined,c);const b=el('div',undefined,'bar',c);el('span',undefined,undefined,b).style.width=pct(v.coverage);}
function showUnified(r){
 const p=$('unified');p.replaceChildren();const h=r.hierarchical;
 const hero=el('article',undefined,'panel hero',p),left=el('div',undefined,undefined,hero),right=el('div',undefined,undefined,hero);
 el('div','LAYER 3 · INTEGRATED RESEARCH INDEX','layerlabel',left);el('h3','综合研究指数',undefined,left);
 el('div',h.reportable?bounds(h.total):'尚未定量','big',left);el('p','量化覆盖 '+pct(h.total.coverage)+' · 功能实验家族 '+h.functional_family_n+'/3',undefined,left);
 el('span',r.unified.evidence,'badge',left);el('h3',r.unified.summary,undefined,right);
 el('p',r.normalization.hgvs_p||r.normalization.hgvs_c,undefined,right);el('p','适用范围：'+readable(r.unified.applicability),undefined,right);
 if(!h.reportable)el('p','该记录不满足量化覆盖与实验家族条件。统一结论保留机制与证据等级，不赋零或满分。',undefined,right);
 else el('p','显示固定权重下兼容的 0–100 区间。缺失项允许 0–100，多来源记录保留范围；不取区间中点作为实测值。',undefined,right);
 if(h.allele_status==='hold_identity_or_splicing')notice(p,'此指数仅描述假设蛋白。身份或剪接门控尚未解决，不能作为该核苷酸等位基因的有效评分。');
 const row=el('div',undefined,'levelrow',p);levelCard(row,'结构 S',h.structure);levelCard(row,'功能缺失 F',h.function);levelCard(row,'显性负性 D',h.dne);
 const fam=el('article',undefined,'panel subpanel',p);el('div','LAYER 2 · FAMILY-BALANCED AGGREGATION','layerlabel',fam);el('h3','每一层怎样汇总',undefined,fam);
 table(fam,['层级 / 实验家族','区间','覆盖','层内规则'],[
 ['结构',h.structure.coverage?bounds(h.structure):'NA',pct(h.structure.coverage),'apo 与锌两项等权；DNA 原值单列'],
 ...h.function.components.map(c=>[c.name,c.interval.coverage?bounds(c.interval):'NA',pct(c.interval.coverage),'各实验家族 1/3；Giacomelli 两个 LOF 条件先等权']),
 ['显性负性',h.dne.coverage?bounds(h.dne):'NA',pct(h.dne.coverage),'仅 Giacomelli WT/nutlin 信号']]);
 el('p',`总指数 = ${num(weights.structure)}·S + ${num(weights.function)}·F + ${num(weights.dne)}·D，再除以权重总和。`,'hint',fam);
 if(h.total.coverage>0)el('p','已覆盖项条件均值范围：'+num(h.total.observed_low)+'–'+num(h.total.observed_high)+'；它不替代包含缺失的总区间。','hint',fam);
 notice(p,'研究指数尚未经 PDO、RNAseq 或临床外部校准。百分位没有统一生物学单位；覆盖率表示量化完整度，不是正确率或证据可信概率。');
 if(r.conflicts.length)r.conflicts.forEach(c=>notice(p,c));
 if(r.manual.length){const m=el('article',undefined,'panel subpanel',p);el('h3','补充实验记录 · 未独立核验',undefined,m);table(m,['维度','指标','值','单位','条件','来源'],r.manual.map(v=>[readable(v.dimension),v.metric,v.value,v.unit,v.condition,v.source]));el('p','缺少预先定义的可比参照，未加入研究指数。','hint',m);}
}
function showSources(r){
 const p=$('sources');p.replaceChildren();p.classList.add('sources');const a=el('article',undefined,'panel',p);el('div','LAYER 1 · SOURCE VALUES AND RANKS','layerlabel',a);el('h3','原始值与统一方向',undefined,a);
 table(a,['指标','原值 / 范围','单位','缺陷方向百分位','参考 N'],db.axis_specs.map((v,i)=>[v.label,r.raw[i],v.unit,r.ranks[i],v.reference_n]));
 const k=el('article',undefined,'panel subpanel',p);el('h3','Kato 八启动子原值',undefined,k);table(k,['启动子','% WT'],PROMOTERS.map((v,i)=>[v,r.promoters[i]]));
 if(r.kotler_records.length){const d=el('details',undefined,undefined,p);el('summary','Kotler 来源记录：'+r.kotler_records.length+' 条',undefined,d);table(d,['accession','RFS'],r.kotler_records);el('p','多记录的核苷酸背景未确定；保持逐记录与观察范围，不赋任意点估计。','hint',d);}
 const classes=el('article',undefined,'panel subpanel',p);el('h3','来源分类 / 计算预测',undefined,classes);
 table(classes,['来源','评估'],[['Kato',r.classes.kato_class],['Giacomelli LOF 策展',r.classes.giacomelli_class_curated],['NCI DNE/LOF',r.classes.dne_lof_class_nci],['Kotler 策展',r.classes.kotler_class_curated],['ClinGen 初步功能代码（未应用）',r.classes.clingen_preliminary_functional_code],['AlphaMissense 分值 / 上游类别',num(r.classes.alphamissense_score)+' / '+num(r.classes.alphamissense_class)]]);
 el('p','AlphaMissense 不参与实验综合指数；ClinGen 初步功能代码不等于完整致病性判定。','hint',classes);
 for(const [name,s] of Object.entries(db.sources)){const d=el('details',undefined,undefined,p);el('summary',name+' · '+s.assessment_kind,undefined,d);el('p',s.assay_context||'历史策展来源','hint',d);ref(d,'原始来源 ↗',s.source.reference);el('small',JSON.stringify(s.source.version)+' · SHA256 '+s.source.sha256,undefined,d);}
}
function showMechanism(r){
 const p=$('mechanism');p.replaceChildren();const s=r.structure,n=r.non_missense;
 const a=el('article',undefined,'panel',p);el('h3','具体替换的结构证据',undefined,a);
 if(s){for(const c of s.categories||[]){el('p',c.label+' · '+c.basis,undefined,a);ref(a,'原始研究 ↗',c.reference);}for(const role of s.residue_roles||[])el('p',role.role+' · 残基层面背景','hint',a);
 table(a,['指标','原始值','单位 / 误差'],(s.measurements||[]).map(m=>[m.metric,m.value,m.unit+(m.error!=null?'; ±'+num(m.error)+' '+m.error_type:'')]));
 if(s.dna_binding?.length){const d=el('details',undefined,undefined,a);el('summary','DNA 结合：10 个识别元件，4°C',undefined,d);table(d,['元件','Kd (µM)','拟合 SE','重复数'],s.dna_binding.map(m=>[m.recognition_element,m.lower_bound!==null?'>'+m.lower_bound:m.value,m.error_se,m.replicate_n]));el('p','>25 为右删失下界，未当作精确 25 或参与雷达 / 总分。','hint',d);}
 if(s.classification_conflict)notice(a,s.classification_conflict);if(s.zinc_classification_caveat)notice(a,s.zinc_classification_caveat);
 }else el('p','本结构来源未覆盖该具体替换或不适用于本变异类型。','hint',a);
 el('p','结构原值来自纯化 DBD；apo 与锌为 10°C，DNA 为 4°C。不能直接外推全长蛋白在体温下的状态。','hint',a);
 const b=el('article',undefined,'panel subpanel',p);el('h3','非 missense / 蛋白底物',undefined,b);
 table(b,['项目','判读'],[['变异类型',readable(n.variant_type)],['序列后果',readable(n.sequence_consequence)],['翻译转录本',r.normalization.translation_transcript],['PTC（变异序列编码位置）',n.ptc_c_position],['PTC 外显子',n.ptc_exon],['NMD 位置预测',n.nmd.predicted===true?'支持，未实测':n.nmd.predicted===false?'不支持，未实测':'未决 / 不适用'],['NMD 状态',readable(n.nmd.status)],['蛋白测量',readable(n.protein_presence)],['重折叠阴性对照',readable(n.refolding_negative_control)]]);
 if(n.domains_if_translated.length)table(b,['若翻译的结构域','原始序列保留'],n.domains_if_translated.map(m=>[readable(m.domain),readable(m.status)]));n.warnings.forEach(v=>notice(b,v));
 if(n.protein_evidence)notice(b,'用户未核验记录：'+n.protein_evidence.assay+' · '+n.protein_evidence.reference+'。需核对灵敏度、多表位和其它等位基因。');
 const c=el('article',undefined,'panel subpanel',p);el('h3','可验证的研究假说',undefined,c);el('p','APR-246/MQ：DNA 接触型不能作为无效硬判据；区分蛋白恢复与红氧作用。','hint',c);ref(c,'MQ 原始研究 ↗','https://www.nature.com/articles/s41419-018-0463-7');
 el('p','ZMC1：锌缺陷提供机制实验线索，但不同突变的救援不一致，不能据此判定临床优于 APR-246。','hint',c);ref(c,'结构与救援原始研究 ↗','https://doi.org/10.7554/eLife.61487');
 el('p','截断或剪接变异：先确认 RNA 产物、蛋白和目标结构；无适合底物只限制蛋白依赖的重折叠假说。','hint',c);ref(c,'TP53 翻译再起始研究 ↗','https://pmc.ncbi.nlm.nih.gov/articles/PMC6614817/');
}
function select(i){current=i;$('selected').value=i;document.querySelectorAll('.variantcard').forEach((v,k)=>v.classList.toggle('active',k===i));showUnified(results[i]);showSources(results[i]);showMechanism(results[i]);}
function render(){
 $('overview').replaceChildren();$('selected').replaceChildren();results.forEach((r,i)=>{const name=r.normalization.hgvs_p||r.normalization.hgvs_c,card=el('article',undefined,'variantcard',$('overview'));card.style.setProperty('--color',colors[i]);card.setAttribute('role','button');card.tabIndex=0;card.setAttribute('aria-label','查看 '+name);el('h3',name,undefined,card);el('strong',r.hierarchical.reportable?bounds(r.hierarchical.total):'尚未定量',undefined,card);el('p','覆盖 '+pct(r.hierarchical.total.coverage)+' · '+readable(r.normalization.kind),undefined,card);el('p',r.unified.summary,undefined,card);card.addEventListener('click',()=>select(i));card.addEventListener('keydown',e=>{if(['Enter',' '].includes(e.key)){e.preventDefault();select(i);}});const o=el('option',name,undefined,$('selected'));o.value=i;});
 $('result-title').textContent=results.length+' 个变异 · 统一评估';$('results').hidden=false;select(Math.min(current,results.length-1));showCharts();
}
function readManual(){const list=[];for(const row of $('measurement-rows').children){const m={};row.querySelectorAll('[data-key]').forEach(i=>m[i.dataset.key]=i.value.trim());if(['metric','value','unit','condition','source'].every(k=>!m[k]))continue;if(!m.value)throw Error('补充测量需填写数值');m.value=Number(m.value);list.push(m);}return list;}
function addMeasurement(){if($('measurement-rows').children.length>=20)return;const row=el('div',undefined,'manualrow',$('measurement-rows')),f=el('div',undefined,'fields',row);for(const [k,t] of [['dimension','维度'],['metric','指标'],['value','原值'],['unit','单位'],['condition','条件'],['source','来源']]){const label=el('label',t,undefined,f),i=el(k==='dimension'?'select':'input',undefined,undefined,label);i.dataset.key=k;if(k==='dimension'){for(const [v,t] of [['structure','结构'],['function','功能']]){const o=el('option',t,undefined,i);o.value=v;}}else if(k==='value'){i.type='number';i.step='any';}}const b=el('button','移除此项',undefined,row);b.type='button';b.addEventListener('click',()=>row.remove());}
function submit(event){event?.preventDefault();$('error').hidden=true;$('results').hidden=true;try{if(!db)throw Error('公开数据尚未校验完成');const q=$('variants').value.split(/[,，\n]+/).map(v=>v.trim()).filter(Boolean);if(!q.length||q.length>6)throw Error('请输入 1–6 个变异');const context={protein_presence:$('protein').value,protein_assay:$('assay').value,protein_evidence_reference:$('protein-source').value,splicing_ruled_out:$('splicing').checked,dimension_inputs:readManual()};if(q.length>1&&(context.protein_presence!=='unmeasured'||context.dimension_inputs.length||context.splicing_ruled_out))throw Error('补充记录与背景请在单变异查询中填写，以免误配其它变异');const next=q.map(v=>interpret(v,db,context,weights));queries=q;lastContext=context;results=next;current=0;render();}catch(e){$('error').textContent=e.message;$('error').hidden=false;}}
function weightsChanged(){try{const v=Object.fromEntries(['structure','function','dne'].map(k=>[k,Number($('w-'+k).value)]));const sum=Object.values(v).reduce((a,b)=>a+b,0);if(!sum)throw Error('至少保留一项非零权重');weights=v;for(const k in v)$('w-'+k+'-out').value=Math.round(v[k]/sum*100)+'%';if(results.length&&!$('results').hidden){results=queries.map(q=>interpret(q,db,lastContext,weights));render();}$('error').hidden=true;}catch(e){$('error').textContent=e.message;$('error').hidden=false;}}
function markdown(){const r=results[current],h=r.hierarchical;return '# TP53 统一研究评估：'+(r.normalization.hgvs_p||r.normalization.hgvs_c)+'\n\n'+r.unified.summary+'\n\n- 综合研究指数：'+(h.reportable?bounds(h.total):'未量化')+'\n- 加权覆盖：'+pct(h.total.coverage)+'\n- 结构：'+levelValue(h.structure)+'\n- 功能缺失：'+levelValue(h.function)+'\n- DNE：'+levelValue(h.dne)+'\n- 权重：'+JSON.stringify(weights)+'\n- 适用：'+readable(r.unified.applicability)+'\n- 软件：'+r.version+'\n- 来源快照：'+r.snapshot_id+'\n\n缺失兼容区间不是置信区间；预设权重尚未外部校准，分值不是致病概率或治疗响应概率。\n\n## 原始指标\n\n|指标|原值|缺陷方向百分位|N|\n|---|---|---|---|\n'+db.axis_specs.map((a,i)=>'|'+[a.label,rank(r.raw[i]),rank(r.ranks[i]),a.reference_n].join('|')+'|').join('\n')+'\n\n## 完整结构与序列机制\n\n```json\n'+JSON.stringify({structure:r.structure,sequence:r.non_missense,manual:r.manual},null,2)+'\n```\n';}
$('query').addEventListener('submit',submit);$('selected').addEventListener('change',()=>select(Number($('selected').value)));$('add-measurement').addEventListener('click',addMeasurement);addMeasurement();
document.querySelectorAll('[data-query]').forEach(b=>b.addEventListener('click',()=>{$('variants').value=b.dataset.query;submit();}));
for(const k of ['structure','function','dne'])$('w-'+k).addEventListener('input',weightsChanged);
document.querySelectorAll('[data-weights]').forEach(b=>b.addEventListener('click',()=>{const v=b.dataset.weights.split(',').map(Number);['structure','function','dne'].forEach((k,i)=>$('w-'+k).value=v[i]);weightsChanged();}));
document.querySelectorAll('[data-tab]').forEach(b=>b.addEventListener('click',()=>{document.querySelectorAll('[data-tab]').forEach(v=>v.setAttribute('aria-selected',v===b));document.querySelectorAll('.tabview').forEach(v=>v.hidden=v.id!==b.dataset.tab);}));
$('download-json').addEventListener('click',()=>save(JSON.stringify(results,null,2),'application/json','TP53_unified_v04.json'));$('download-svg').addEventListener('click',()=>save(new XMLSerializer().serializeToString($('radar').querySelector('svg')),'image/svg+xml','TP53_radar_v04.svg'));$('download-report').addEventListener('click',()=>save(markdown(),'text/markdown;charset=utf-8','TP53_unified_v04.md'));
async function load(){try{
 $('submit').disabled=true;const [d,m]=await Promise.all([fetch('./data/public.json.gz'),fetch('./data/manifest.json')]);if(!d.ok||!m.ok)throw Error('数据文件读取失败');const buffer=await d.arrayBuffer(),magic=new Uint8Array(buffer),manifest=await m.json();
 const text=magic[0]===31&&magic[1]===139?await new Response(new Blob([buffer]).stream().pipeThrough(new DecompressionStream('gzip'))).text():new TextDecoder().decode(buffer);
 const hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(text)))).map(x=>x.toString(16).padStart(2,'0')).join('');if(hash!==manifest.public_json_sha256)throw Error('数据校验和不一致，停止解读');
 db=JSON.parse(text);const profiles=Object.values(db.variants).map(v=>scoreProfile(v.ranks)),eligible=profiles.filter(v=>v.reportable).length,complete=profiles.filter(v=>v.total.coverage>=1-1e-10).length;
 $('coverage-summary').textContent='默认权重下：'+eligible.toLocaleString()+'/'+profiles.length.toLocaleString()+' 个蛋白替换满足总指数显示条件；'+complete+' 个全部分量有测量。其它变异保留分层证据和统一定性结论。';
 $('health').textContent=Object.keys(db.variants).length.toLocaleString()+' 种替换 · 快照已校验';$('build').textContent='SOURCE '+db.snapshot_id.slice(0,12)+' · SCORE research-index-0.4.0';$('submit').disabled=false;submit();
 }catch(e){$('health').textContent='数据校验失败';$('error').textContent=e.message;$('error').hidden=false;}}
load();
