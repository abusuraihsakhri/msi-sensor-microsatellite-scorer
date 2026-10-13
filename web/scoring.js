"use strict";
function scoreMarkers(markers){
 const names=["BAT-25","BAT-26","D2S123","D5S346","D17S250"];
 if(!markers||typeof markers!=="object")throw Error("Five-marker panel is required");
 const keys=Object.keys(markers);
 if(keys.length!==names.length||names.some(x=>!Object.prototype.hasOwnProperty.call(markers,x)||typeof markers[x]!=="boolean"))throw Error("Provide exactly five stable/unstable marker results");
 const unstable=names.filter(x=>markers[x]);
 return {msi_status:unstable.length>=2?"MSI-H":unstable.length===1?"MSI-L":"MSS",unstable_count:unstable.length,total_markers:5,unstable_markers:unstable,percent_unstable:unstable.length*20};
}

function interpretIHC(proteins, braf=null, methylation=null){
 const names=["MLH1","MSH2","MSH6","PMS2"];
 if(!proteins||names.some(x=>typeof proteins[x]!=="boolean")||Object.keys(proteins).length!==4)throw Error("Provide all four retained/lost IHC values");
 if((braf!==null&&typeof braf!=="boolean")||(methylation!==null&&typeof methylation!=="boolean"))throw Error("Invalid optional IHC molecular result");
 const lost=names.filter(x=>!proteins[x]);
 let etiology="No MMR protein loss identified",lynch=false;
 if(lost.length){
  lynch=null;etiology="Unusual loss pattern; specialist evaluation recommended";
  if(!proteins.MLH1&&!proteins.PMS2){
   if(methylation===true){etiology="MLH1 methylation detected; sporadic origin favored";lynch=false;}
   else if(braf===true){etiology="BRAF V600E detected; sporadic origin favored in colorectal tumors";lynch=false;}
   else if(braf===false&&methylation===false){etiology="Lynch syndrome evaluation warranted";lynch=true;}
   else{etiology="MLH1/PMS2 loss; further molecular testing indicated";}
  }else if(!proteins.MSH2&&!proteins.MSH6){etiology="MSH2/MSH6 loss; Lynch syndrome evaluation warranted";lynch=true;}
  else if(!proteins.MSH6&&proteins.MSH2){etiology="Isolated MSH6 loss; Lynch syndrome evaluation warranted";lynch=true;}
  else if(!proteins.PMS2&&proteins.MLH1){etiology="Isolated PMS2 loss; Lynch syndrome evaluation warranted";lynch=true;}
 }
 return {mmr_status:lost.length?"dMMR":"pMMR",deficient_proteins:lost,retained_proteins:names.filter(x=>proteins[x]),etiology,suspected_lynch_syndrome:lynch};
}

function combined(markers=null,proteins=null,braf=null,methylation=null){
 if(!markers&&!proteins)throw Error("Provide PCR markers or IHC results");
 const pcr=markers?scoreMarkers(markers):null,ihc=proteins?interpretIHC(proteins,braf,methylation):null;
 let concordance=null,final_status=null;
 if(pcr&&ihc){
  concordance=(pcr.msi_status==="MSI-H")===(ihc.mmr_status==="dMMR")?"Concordant":"Discordant";
  final_status=concordance==="Discordant"?"Discordant (PCR/IHC)":pcr.msi_status==="MSI-H"?"MSI-H/dMMR":pcr.msi_status==="MSI-L"?"MSI-L/pMMR":"MSS/pMMR";
 }else final_status=pcr?pcr.msi_status:ihc.mmr_status;
 return {msi_pcr_result:pcr,mmr_ihc_result:ihc,concordance,final_status};
}

function parseCSV(text){
 const rows=[];let row=[],field="",quoted=false,afterQuote=false;
 for(let i=0;i<text.length;i++){
  const ch=text[i];
  if(quoted){if(ch==='"'){if(text[i+1]==='"'){field+='"';i++;}else{quoted=false;afterQuote=true;}}else field+=ch;continue;}
  if(afterQuote){if(ch!==","&&ch!=="\n"&&ch!=="\r")throw Error("Unexpected character after CSV quote");afterQuote=false;}
  if(ch==='"'){if(field)throw Error("Unexpected CSV quote");quoted=true;}
  else if(ch===","){row.push(field);field="";}
  else if(ch==="\n"||ch==="\r"){if(ch==="\r"&&text[i+1]==="\n")i++;row.push(field);rows.push(row);row=[];field="";}
  else field+=ch;
 }
 if(quoted)throw Error("Unterminated CSV quoted field");
 if(field!==""||row.length){row.push(field);rows.push(row);}
 return rows;
}

function batchCSV(text){
 if(text.length>2_000_000)throw Error("CSV exceeds 2 MB browser processing limit");
 const rows=parseCSV(text.replace(/^\uFEFF/,""));
 if(!rows.length)throw Error("CSV has no header");
 const head=rows[0],names=["BAT-25","BAT-26","D2S123","D5S346","D17S250"],mmr=["MLH1","MSH2","MSH6","PMS2"];
 if(new Set(head).size!==head.length||names.some(x=>!head.includes(x)))throw Error("CSV needs unique headers and all five Bethesda markers");
 const hasIHC=mmr.some(x=>head.includes(x));
 if(hasIHC&&mmr.some(x=>!head.includes(x)))throw Error("IHC requires MLH1, MSH2, MSH6 and PMS2 columns");
 const extra=["msi_status","unstable_count","mmr_status","final_status"];
 if(extra.some(x=>head.includes(x)))throw Error("CSV contains reserved output columns");
 if(rows.length>5001)throw Error("Limit: 5,000 records per CSV");
 const parseBool=(value,kind,label,line)=>{
  const t=(value??"").trim().toLowerCase();
  const yes=kind==="ihc"?["true","1","yes","retained","positive","present","detected","methylated"]:["true","1","yes","unstable"];
  const no=kind==="ihc"?["false","0","no","lost","negative","absent","not_detected","unmethylated"]:["false","0","no","stable"];
  if(yes.includes(t))return true;
  if(no.includes(t))return false;
  throw Error("Row "+line+": invalid "+label+" value");
 };
 const escapeValue=v=>{
  let s=String(v??"");
  if(/^\s*[=+\-@]/.test(s))s="'"+s;
  return /[",\r\n]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s;
 };
 const out=[head.concat(extra).map(escapeValue).join(",")];
 const index=name=>head.indexOf(name);
 for(let n=1;n<rows.length;n++){
  const row=rows[n];if(row.length===1&&row[0]==="")continue;
  if(row.length!==head.length)throw Error("Row "+(n+1)+": wrong number of columns");
  const markers=Object.fromEntries(names.map(x=>[x,parseBool(row[index(x)],"pcr",x,n+1)]));
  const ihc=hasIHC?Object.fromEntries(mmr.map(x=>[x,parseBool(row[index(x)],"ihc",x,n+1)])):null;
  const opt=col=>{const ix=index(col);return ix<0||!row[ix].trim()?null:parseBool(row[ix],"ihc",col,n+1);};
  const result=combined(markers,ihc,hasIHC?opt("BRAF_V600E"):null,hasIHC?opt("MLH1_PROMOTER_METHYLATION"):null);
  out.push(row.concat([result.msi_pcr_result.msi_status,result.msi_pcr_result.unstable_count,result.mmr_ihc_result?.mmr_status??"N/A",result.final_status]).map(escapeValue).join(","));
 }
 return {csv:out.join("\r\n")+"\r\n",count:out.length-1};
}

const MSIScorer={scoreMarkers,interpretIHC,combined,parseCSV,batchCSV};
if(typeof module!=="undefined"&&module.exports)module.exports=MSIScorer;
if(typeof window!=="undefined")window.MSIScorer=MSIScorer;
