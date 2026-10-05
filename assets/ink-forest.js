'use strict';
'use strict';
// =====================================================================================================================================
// Ink Forest: an interactive Lya tomography game. https://jibancat.github.io/ink-forest/  (code: MIT licence; data: see DATA.md)
// The scientific object is a 2D tomography sheet through one simulated (PRIYA) plane of the intergalactic medium at z = 3:
// x = depth along the line of sight (drawn downward), y = transverse position. Each observed vertical line is that sightline's own
// noisy synthetic Lya forest spectrum; the field between lines is the exact Gaussian conditional mean of the smoothed absorption
// (1.5 Mpc/h) given the observed lines only. Nothing in this page knows the true field between the lines. See SCIENCE.md.
// =====================================================================================================================================
/* the public build takes no URL flags */const Q_=new URLSearchParams(''),MODE=Q_.get('mode')||'play';let REPV=(Q_.get('rep')||'A').toUpperCase();
const D=window.SHEETDATA,M=D.meta,N=M.N,NQ=M.Q+1,WINB=M.win_bins,ALL=M.rows,LAG=M.lag_rows,NPXW=M.npx_win;
function b64(s,T){const b=atob(s),u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);return new T(u.buffer);}
const OBS=b64(D.obs,Int16Array),CDD=b64(D.cdd,Float32Array),CTD=b64(D.ctd,Float32Array);
const NPXF=M.npx_full,FLUXF=b64(D.fluxfull,Uint8Array),OBSF=b64(D.obsfull,Int16Array),CDDF=b64(D.cddfull,Float32Array),NQF=M.q_full+1;
// ---- the LOS window: one of three declared 40 Mpc/h tiles of the same sightlines (0-40, 40-80, 80-120 Mpc/h); the game shows 40-80.
const X0Q=Q_.get('x0'),X0M=X0Q!=null&&[0,40,80].includes(+X0Q)?+X0Q:40,X0B=X0M*4,X0P=Math.round(X0M/M.px_Mpc_h);
// ---- the survey: NAV candidate lines every 6 rows (1.5 Mpc/h) centred on the survey region; the sheet spans half a spacing beyond them
const NAV=Math.max(4,Math.min(25,+(Q_.get('n')||25))),STEP=M.step_rows,R0=335-STEP*(NAV>>1),SURV=Array.from({length:NAV},(_,k)=>R0+STEP*k);
const LO=SURV[0]-(STEP>>1),NR=STEP*NAV;   // sheet rows LO .. LO+NR-1 (0.25 Mpc/h each)
const ROWIDX=r=>ALL.indexOf(r);
// ---- the estimator in the plane (exact per-LOS-mode solve, modes q <= Q): Gaussian conditional mean and variance (see SCIENCE.md)
const DRE=new Float64Array(ALL.length*NQ),DIM=new Float64Array(ALL.length*NQ);
for(let k=0;k<ALL.length;k++)for(let q=0;q<NQ;q++){let s=0,t=0;for(let x=0;x<N;x++){const v=OBS[k*N+x]*M.obs_scale,a=2*Math.PI*q*x/N;s+=v*Math.cos(a);t-=v*Math.sin(a);}DRE[k*NQ+q]=s;DIM[k*NQ+q]=t;}
const COS=new Float64Array(NQ*WINB),SIN=new Float64Array(NQ*WINB);for(let q=0;q<NQ;q++)for(let x=0;x<WINB;x++){COS[q*WINB+x]=Math.cos(2*Math.PI*q*(x+X0B)/N);SIN[q*WINB+x]=Math.sin(2*Math.PI*q*(x+X0B)/N);}
let lastSolveMs=0;
// ---- the damped-absorber rule (unchanged from earlier Ink Forest versions), from the line's OWN observed spectrum:
// pixel pairs -> 694 samples of 20 km/s; A = clamp(1 - mean over +-2 samples); a damped core is a run of A > 0.93 of >= 40 samples; the
// mask widens the core along the line's own absorption while A > 0.3. Mask bins on the estimator's 0.25 Mpc/h grid.
const ABS={};
function absorberOf(r){if(r in ABS)return ABS[r];const k=ROWIDX(r),NS=694,Fn=new Float64Array(NS),A=new Float64Array(NS);if(k<0)return ABS[r]=null;
  for(let i=0;i<NS;i++)Fn[i]=0.5*(OBSF[k*NPXF+2*i]+OBSF[k*NPXF+2*i+1])/20000;
  for(let i=0;i<NS;i++){let s_=0,c_=0;for(let j=Math.max(0,i-2);j<Math.min(NS,i+3);j++){s_+=Fn[j];c_++;}A[i]=Math.min(1,Math.max(0,1-s_/c_));}
  let dla=null;for(let i=0;i<NS;){if(A[i]>0.93){let j=i;while(j<NS&&A[j]>0.93)j++;if(j-i>=40){let a=i,b=j-1;while(a>0&&A[a-1]>0.3)a--;while(b<NS-1&&A[b+1]>0.3)b++;dla={i0:a,i1:b,c0:i,c1:j-1};}i=j;}else i++;}
  if(dla){const SM=2*M.px_Mpc_h;dla.b0=Math.floor(dla.i0*SM/0.25+1e-9);dla.b1=Math.ceil((dla.i1+1)*SM/0.25-1e-9)-1;dla.Mpc=[dla.i0*SM,(dla.i1+1)*SM];}
  return ABS[r]=dla;}
// ---- the estimator with set-aside lines masked (exact): K'^-1 = K^-1 - K^-1 P (P^T K^-1 P)^-1 P^T K^-1 over ALL LOS modes
// (validated against a dense pixel-space solve to 1e-14). The masked bins' data cannot influence the result; support becomes c(t, x).
let DREF=null,DIMF=null,COSN=null;
function fullModes(){if(DREF)return;COSN=new Float64Array(N);for(let j=0;j<N;j++)COSN[j]=Math.cos(2*Math.PI*j/N);DREF=new Float64Array(ALL.length*NQF);DIMF=new Float64Array(ALL.length*NQF);
  for(let k=0;k<ALL.length;k++)for(let q=0;q<NQF;q++){let a=0,b=0;for(let x=0;x<N;x++){const v=OBS[k*N+x]*M.obs_scale,j=(q*x)%N;a+=v*COSN[j];b-=v*COSN[(j-N/4+N)%N];}DREF[k*NQF+q]=a;DIMF[k*NQF+q]=b;}}
const cosN=j=>COSN[((j%N)+N)%N],sinN=j=>COSN[(((j-N/4)%N)+N)%N];
// The exact masked solve runs as a resumable job of small, deterministic units (a generator), so no frame carries it: setup (all-mode
// inverses, the masked-pixel system), then the field (mean) for every row, then the exact support row by row. Same arithmetic as a single pass.
function* maskedSteps(rows,mrows,out){
  if(!DREF){COSN=new Float64Array(N);for(let j=0;j<N;j++)COSN[j]=Math.cos(2*Math.PI*j/N);DREF=new Float64Array(ALL.length*NQF);DIMF=new Float64Array(ALL.length*NQF);
    for(let k=0;k<ALL.length;k++){for(let q=0;q<NQF;q++){let a=0,b=0;for(let x=0;x<N;x++){const v=OBS[k*N+x]*M.obs_scale,j=(q*x)%N;a+=v*COSN[j];b-=v*COSN[(j-N/4+N)%N];}DREF[k*NQF+q]=a;DIMF[k*NQF+q]=b;}yield 1;}}
  const n=rows.length,ki=rows.map(ROWIDX),{s,c}=out,Ki=new Float64Array(NQF*n*n),WR=new Float64Array(NQF*n),WI=new Float64Array(NQF*n),K=new Float64Array(n*n),Lc=new Float64Array(n*n);
  for(let q=0;q<NQF;q++){for(let i=0;i<n;i++)for(let j=0;j<n;j++)K[i*n+j]=CDDF[Math.abs(rows[i]-rows[j])*NQF+q]+(i===j?M.noise_var:0);
    for(let i=0;i<n;i++)for(let j=0;j<=i;j++){let v=K[i*n+j];for(let k=0;k<j;k++)v-=Lc[i*n+k]*Lc[j*n+k];Lc[i*n+j]=i===j?Math.sqrt(Math.max(v,1e-30)):v/Lc[j*n+j];}
    const o=q*n*n,y=new Float64Array(n),z=new Float64Array(n);
    for(let col=0;col<n;col++){for(let i=0;i<n;i++){let v=i===col?1:0;for(let k=0;k<i;k++)v-=Lc[i*n+k]*y[k];y[i]=v/Lc[i*n+i];}
      for(let i=n-1;i>=0;i--){let v=y[i];for(let k=i+1;k<n;k++)v-=Lc[k*n+i]*z[k];z[i]=v/Lc[i*n+i];}for(let i=0;i<n;i++)Ki[o+i*n+col]=z[i];}
    for(let i=0;i<n;i++){let a=0,b=0;for(let j=0;j<n;j++){a+=Ki[o+i*n+j]*DREF[ki[j]*NQF+q];b+=Ki[o+i*n+j]*DIMF[ki[j]*NQF+q];}WR[q*n+i]=a;WI[q*n+i]=b;}
    if(q%60===59)yield 1;}
  // masked pixels (line index, LOS bin)
  const PK=[],PX=[];for(const r of mrows){const d=absorberOf(r),k=rows.indexOf(r);for(let b=d.b0;b<=d.b1;b++){PK.push(k);PX.push(b);}}const m=PK.length;
  const mk=[...new Set(PK)],gab={};const wq=q=>(q===0||q===N/2)?1:2;
  for(const a of mk)for(const b of mk){const g=new Float64Array(N);for(let dl=0;dl<N;dl++){let v=0;for(let q=0;q<NQF;q++)v+=wq(q)*Ki[q*n*n+a*n+b]*cosN(q*dl);g[dl]=v/N;if(dl%120===119)yield 1;}gab[a+','+b]=g;}
  const G=new Float64Array(m*m),LG=new Float64Array(m*m);for(let i=0;i<m;i++)for(let j=0;j<m;j++)G[i*m+j]=gab[PK[i]+','+PK[j]][((PX[i]-PX[j])%N+N)%N];
  for(let i=0;i<m;i++){for(let j=0;j<=i;j++){let v=G[i*m+j];for(let k=0;k<j;k++)v-=LG[i*m+k]*LG[j*m+k];LG[i*m+j]=i===j?Math.sqrt(Math.max(v,1e-30)):v/LG[j*m+j];}if(i%50===49)yield 1;}
  const fwd=(v)=>{for(let i=0;i<m;i++){let a=v[i];const o=i*m;for(let k=0;k<i;k++)a-=LG[o+k]*v[k];v[i]=a/LG[o+i];}return v;};
  const bwd=(v)=>{for(let i=m-1;i>=0;i--){let a=v[i];for(let k=i+1;k<m;k++)a-=LG[k*m+i]*v[k];v[i]=a/LG[i*m+i];}return v;};
  const u=new Float64Array(m);for(let i=0;i<m;i++){const k=PK[i];let v=0;for(let q=0;q<NQF;q++)v+=wq(q)*(WR[q*n+k]*cosN(q*PX[i])-WI[q*n+k]*sinN(q*PX[i]));u[i]=v/N;}bwd(fwd(u));
  for(let q=0;q<NQ;q++){const o=q*n*n;for(const k of mk){let ur=0,ui=0;for(let i=0;i<m;i++)if(PK[i]===k){ur+=u[i]*cosN(q*PX[i]);ui-=u[i]*sinN(q*PX[i]);}
      for(let l=0;l<n;l++){WR[q*n+l]-=Ki[o+l*n+k]*ur;WI[q*n+l]-=Ki[o+l*n+k]*ui;}}}
  yield 1;
  const Mr=new Float64Array(NQ),Mi=new Float64Array(NQ),X=new Float64Array(NQ*n),h=new Float64Array(mk.length*N),V=new Float64Array(m),REDT=new Float64Array(NR);
  // the field (posterior mean) for every row, with red(t)
  for(let t=0;t<NR;t++){const row=LO+t;let red=0;
    for(let q=0;q<NQ;q++){let r_=0,mi=0;for(let l=0;l<n;l++){const v=CTD[Math.abs(row-rows[l])*NQ+q];X[q*n+l]=v;r_+=v*WR[q*n+l];mi+=v*WI[q*n+l];}Mr[q]=r_;Mi[q]=mi;
      const o=q*n*n;let tt=0;for(let a=0;a<n;a++){let v=0;for(let b=0;b<n;b++)v+=Ki[o+a*n+b]*X[q*n+b];tt+=X[q*n+a]*v;}red+=(q?2:1)*tt;}
    for(let x=0;x<WINB;x++){let v=Mr[0];for(let q=1;q<NQ;q++)v+=2*(Mr[q]*COS[q*WINB+x]-Mi[q]*SIN[q*WINB+x]);s[t*WINB+x]=v/N/M.sigmaT;}
    REDT[t]=red;if(t%25===24)yield 1;}
  out.meanReady=true;yield 1;
  // support: red(t) minus the information the masked bins would have carried, v^T G^-1 v with v_i = h_t(k_i; x_i - x)
  for(let t=0;t<NR;t++){const row=LO+t,red=REDT[t];for(let q=0;q<NQ;q++)for(let l=0;l<n;l++)X[q*n+l]=CTD[Math.abs(row-rows[l])*NQ+q];
    mk.forEach((k,ii)=>{const zq=new Float64Array(NQ);for(let q=0;q<NQ;q++){const o=q*n*n;let v=0;for(let b=0;b<n;b++)v+=Ki[o+k*n+b]*X[q*n+b];zq[q]=v;}
      for(let dl=0;dl<N;dl++){let v=zq[0];for(let q=1;q<NQ;q++)v+=2*zq[q]*cosN(q*dl);h[ii*N+dl]=v/N;}});
    for(let x=0;x<WINB;x++){const xx=X0B+x;for(let i=0;i<m;i++)V[i]=h[mk.indexOf(PK[i])*N+((PX[i]-xx)%N+N)%N];fwd(V);let cr=0;for(let i=0;i<m;i++)cr+=V[i]*V[i];
      c[t*WINB+x]=Math.min(1,Math.max(0,(red/N-cr)/M.varT));}
    yield 1;}}
function posteriorMasked(rows,mrows,lazy){const out={s:new Float32Array(NR*WINB),c:new Float32Array(NR*WINB),masked:mrows.slice(),meanReady:false};
  const it=maskedSteps(rows,mrows,out),job={done:false,units:0,ms:0,meanMs:null,
    // advance up to k units, or until msBudget is spent (live play); count-only in recordings so captures stay deterministic
    run(k,msBudget){const t0=performance.now();for(let i=0;i<k;i++){const r=it.next();this.units++;if(r.done){this.done=true;break;}if(out.meanReady&&this.meanMs==null)this.meanMs=this.ms+performance.now()-t0;if(msBudget&&performance.now()-t0>msBudget)break;}
      this.ms+=performance.now()-t0;return this.done;}};
  if(!lazy)job.run(1e9);
  out.job=lazy&&!job.done?job:null;out.jobInfo=job;return out;}
function posterior(rows,aside,lazy){const mrows=(rows||[]).filter(r=>aside&&(aside.has?aside.has(r):aside.includes(r))&&absorberOf(r));
  if(mrows.length){const t0=performance.now(),f=posteriorMasked(rows,mrows,lazy);lastSolveMs=performance.now()-t0;return f;}return posteriorPlain(rows);}
function posteriorPlain(rows){
  const t0=performance.now(),n=rows.length,s=new Float32Array(NR*WINB),c=new Float32Array(NR);if(!n){lastSolveMs=0;return {s,c};}
  const L=new Float64Array(NQ*n*n),WR=new Float64Array(NQ*n),WI=new Float64Array(NQ*n),K=new Float64Array(n*n),ki=rows.map(ROWIDX);
  for(let q=0;q<NQ;q++){
    for(let i=0;i<n;i++)for(let j=0;j<n;j++)K[i*n+j]=CDD[Math.abs(rows[i]-rows[j])*NQ+q]+(i===j?M.noise_var:0);
    const o=q*n*n;
    for(let i=0;i<n;i++)for(let j=0;j<=i;j++){let v=K[i*n+j];for(let k=0;k<j;k++)v-=L[o+i*n+k]*L[o+j*n+k];L[o+i*n+j]=i===j?Math.sqrt(Math.max(v,1e-30)):v/L[o+j*n+j];}
    for(const [src,dst] of [[DRE,WR],[DIM,WI]]){const y=new Float64Array(n);
      for(let i=0;i<n;i++){let v=src[ki[i]*NQ+q];for(let k=0;k<i;k++)v-=L[o+i*n+k]*y[k];y[i]=v/L[o+i*n+i];}
      for(let i=n-1;i>=0;i--){let v=y[i];for(let k=i+1;k<n;k++)v-=L[o+k*n+i]*dst[q*n+k];dst[q*n+i]=v/L[o+i*n+i];}}}
  const Mr=new Float64Array(NQ),Mi=new Float64Array(NQ),kv=new Float64Array(n),y=new Float64Array(n);
  for(let t=0;t<NR;t++){const row=LO+t;let red=0;
    for(let q=0;q<NQ;q++){let r=0,m=0;const o=q*n*n;
      for(let l=0;l<n;l++){const v=CTD[Math.abs(row-rows[l])*NQ+q];kv[l]=v;r+=v*WR[q*n+l];m+=v*WI[q*n+l];}
      Mr[q]=r;Mi[q]=m;let tt=0;for(let a=0;a<n;a++){let v=kv[a];for(let k=0;k<a;k++)v-=L[o+a*n+k]*y[k];y[a]=v/L[o+a*n+a];tt+=y[a]*y[a];}red+=(q?2:1)*tt;}
    c[t]=Math.min(1,Math.max(0,red/N/M.varT));
    for(let x=0;x<WINB;x++){let v=Mr[0];for(let q=1;q<NQ;q++)v+=2*(Mr[q]*COS[q*WINB+x]-Mi[q]*SIN[q*WINB+x]);s[t*WINB+x]=v/N/M.sigmaT;}}
  lastSolveMs=performance.now()-t0;return {s,c};}
// ---------------------------------------------------------------- GL
const cv=document.getElementById('c'),gl=cv.getContext('webgl2',{alpha:false,antialias:true,preserveDrawingBuffer:MODE==='rec'});
gl.getExtension('EXT_color_buffer_float');
function sh(type,src){const s=gl.createShader(type);gl.shaderSource(s,src);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS)){console.error(src.split('\n').map((l,i)=>(i+1)+': '+l).join('\n'));throw new Error(gl.getShaderInfoLog(s));}return s;}
function prog(vs,fs){const p=gl.createProgram();gl.attachShader(p,sh(gl.VERTEX_SHADER,vs));gl.attachShader(p,sh(gl.FRAGMENT_SHADER,fs));gl.linkProgram(p);if(!gl.getProgramParameter(p,gl.LINK_STATUS))throw new Error(gl.getProgramInfoLog(p));
  const u={};const n=gl.getProgramParameter(p,gl.ACTIVE_UNIFORMS);for(let i=0;i<n;i++){const a=gl.getActiveUniform(p,i);u[a.name.replace('[0]','')]=gl.getUniformLocation(p,a.name);}return {p,u};}
function buf(data,usage=gl.STATIC_DRAW){const b=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.bufferData(gl.ARRAY_BUFFER,data,usage);return b;}
function attrib(pp,name,b,size){const l=gl.getAttribLocation(pp.p,name);if(l<0)return;gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.enableVertexAttribArray(l);gl.vertexAttribPointer(l,size,gl.FLOAT,false,0,0);}
// ---- world geometry (1 unit = 10 Mpc/h): the sheet hangs in the plane z = 0; x across (transverse y), y up; LOS runs down from the top
const SW=NR*0.025,SH=WINB*0.025,SX0=-SW/2,STOP=SH/2;
const rowX=r=>SX0+(r-LO+0.5)*0.025,depthY=b=>STOP-(b+0.5)*0.025;
// ---- field textures (s, c) per sheet cell: width NR (transverse), height WINB (LOS); two for the EMERGE morph
function fieldTex(){const t=gl.createTexture();gl.bindTexture(gl.TEXTURE_2D,t);gl.texImage2D(gl.TEXTURE_2D,0,gl.RG16F,NR,WINB,0,gl.RG,gl.FLOAT,new Float32Array(NR*WINB*2));
  for(const [k,v] of [[gl.TEXTURE_MIN_FILTER,gl.LINEAR],[gl.TEXTURE_MAG_FILTER,gl.LINEAR],[gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE],[gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE]])gl.texParameteri(gl.TEXTURE_2D,k,v);return t;}
const TXF=[fieldTex(),fieldTex()];
const cAt=(f,t,b)=>f.c.length===NR?f.c[t]:f.c[t*WINB+b];   // support is c(t) without a mask, c(t, x) with one
function uploadField(t,f){const d=new Float32Array(NR*WINB*2);for(let r=0;r<NR;r++)for(let x=0;x<WINB;x++){const o=(x*NR+r)*2;d[o]=f.s[r*WINB+x];d[o+1]=cAt(f,r,x);}
  gl.bindTexture(gl.TEXTURE_2D,t);gl.texSubImage2D(gl.TEXTURE_2D,0,0,0,NR,WINB,gl.RG,gl.FLOAT,d);}
// ---- the sheet material (player render). Signed: absorbing excess -> graphite ink; transmissive excess -> lifted, open paper.
// Support c changes MATERIAL CERTAINTY (edge definition, grain, crisp iso-contours), not the presence of the field. Nothing hidden.
const SHEET_VS=`#version 300 es
in vec2 aQ;uniform mat4 uVP;uniform float uFace;uniform vec4 uRect;out vec2 vUV;out vec3 vW;
void main(){vUV=aQ*0.5+0.5;vec3 w=vec3(uRect.x+vUV.x*uRect.z,uRect.y-(1.0-vUV.y)*uRect.w,0.0);vW=w;
 gl_Position=uFace>0.5?vec4(aQ.x,aQ.y,0.0,1.0):uVP*vec4(w,1.0);}`;
const pSheet=prog(SHEET_VS,`#version 300 es
precision highp float;in vec2 vUV;in vec3 vW;uniform sampler2D uF0,uF1;uniform float uMix,uFront,uFrontR,uFrontU,uFrontW,uAlpha,uSwatch,uRep,uTime,uCut,uCutY,uNoField,uLayer;
uniform vec2 uRes;out vec4 o;
float h(vec2 p){return fract(sin(dot(p,vec2(12.9898,78.233)))*43758.5453);}
float vn(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);return mix(mix(h(i),h(i+vec2(1,0)),f.x),mix(h(i+vec2(0,1)),h(i+vec2(1,1)),f.x),f.y);}
vec2 field(vec2 uv){vec2 a=texture(uF0,vec2(uv.x,1.0-uv.y)).rg,b=texture(uF1,vec2(uv.x,1.0-uv.y)).rg;
 float m=uMix;if(uFront>0.5){float d=abs(uv.x-uFrontU)*${(NR*0.25).toFixed(2)};m=smoothstep(d-uFrontW,d,uFrontR);}return mix(a,b,m);}
void main(){vec2 f=uSwatch>0.5?vec2(0.0,1.0):field(vUV);float s=f.r,c=f.g;if(uNoField>0.5){s=0.0;c=1.0;}
 vec2 px=vUV*vec2(${NR}.0,${WINB}.0);
 float fib=vn(px*vec2(1.6,0.45))*0.55+vn(px*3.7)*0.45;
 // grain at a constant ~2 screen px at any distance (two octaves crossfaded, so it neither swims nor turns into coarse static close up)
 float kk=max(max(fwidth(px.x),fwidth(px.y)),1e-4),lv=log2(kk*2.0),l0=floor(lv),fr=lv-l0;
 float tooth=mix(vn(px/exp2(l0)+17.0),vn(px/exp2(l0+1.0)+29.0),fr),tooth2=mix(vn(px/exp2(l0-0.5)+3.0),vn(px/exp2(l0+0.5)+41.0),fr);
 // toned paper (a mid value) so both signs have room: absorbing excess is graphite ink, transmissive excess lifts the paper to clean white
 vec3 PAPER=vec3(0.792,0.782,0.755)*(0.975+0.05*fib),LIFT=vec3(0.985,0.984,0.975),INK=vec3(0.085,0.085,0.10);
 // tone: s in prior sigma -> amount of material, near-linear; the tails are equalised (full ink at s = +1.7, full lift at s = -1.2) because the
 // inferred field is skewed (deep walls, shallow voids). Amplitude is NOT scaled by support. Wash ink stops short of the forest's black.
 float aI=clamp((s-0.10)/1.90,0.0,1.0),aL=clamp((-s-0.08)/1.12,0.0,1.0);
 // material certainty from support c: where c is low the wash is wet (feathered, grainy, soft-edged); where c is high it has dried
 // (smooth, with the pooled rim of dried ink at the strong level, and a clean cut edge round lifted paper). Same mean density either way.
 float wet=1.0-smoothstep(0.30,0.80,c),grain=(tooth-0.5)*0.55+(tooth2-0.5)*0.75;
 aI=clamp(aI*(1.0+wet*grain),0.0,1.0);aL=clamp(aL*(1.0-wet*grain*1.1),0.0,1.0);
 vec3 col=mix(PAPER,LIFT,aL*0.97);col=mix(col,INK,aI*0.80);
 float fw=max(fwidth(s),1e-4),dry=1.0-wet;
 float rimI=smoothstep(0.0,1.0,(s-0.75)/fw)*(1.0-smoothstep(1.0,3.5,(s-0.75)/fw)),rimL=smoothstep(0.0,1.0,(-0.75-s)/fw)*(1.0-smoothstep(1.0,3.0,(-0.75-s)/fw));
 col=mix(col,INK,rimI*0.30*dry);col=mix(col,vec3(1.0),rimL*0.55*dry);
 // the unsurveyed sheet keeps a rough tooth: nobody has measured there (subtle; mean-neutral)
 col*=1.0+wet*(1.0-aI)*(1.0-aL)*0.035*(tooth2-0.5);
 float a=uAlpha;
 if(uLayer>0.5&&uLayer<1.5)col=mix(col,PAPER,0.80);   // debrief 'measured': the inference steps back
 if(uLayer>2.5){float un=1.0-step(0.30,c),hatch=step(0.55,fract((px.x+px.y)*0.25));col=mix(col,mix(col,vec3(0.30,0.30,0.32),0.40*hatch),un);col=mix(col,PAPER,0.45*(1.0-un));}   // 'unresolved'
 if(uRep>1.5&&uCut>0.5){float dy=abs(vW.y-uCutY);col=mix(col,vec3(0.25,0.25,0.27),(1.0-smoothstep(0.0,0.004,dy))*0.6);}
 if(uRep>0.5){a=uAlpha*(1.0-0.55*aL)*(0.84+0.16*aI);}    // curtain: open structure is literally more see-through
 o=vec4(col,a);}`);
const qBuf=buf(new Float32Array([-1,-1,1,-1,-1,1,1,1]));
// ---------------------------------------------------------------- camera
const V3={sub:(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],cross:(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dot:(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],norm:a=>{const l=Math.hypot(...a)||1;return a.map(x=>x/l);},add:(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],mul:(a,k)=>a.map(x=>x*k)};
let VP,EYE,BASIS,CURFOV=0.8;const CAM={pose:{eye:[1.6,0.6,4.2],tgt:[0,0,0],fov:0.8}};
function camUpdate(){const asp=cv.width/cv.height,p=CAM.pose;EYE=p.eye;CURFOV=p.fov;const z=V3.norm(V3.sub(EYE,p.tgt)),x=V3.norm(V3.cross([0,1,0],z)),y=V3.cross(z,x);BASIS={x,y,z};
  const view=[x[0],y[0],z[0],0,x[1],y[1],z[1],0,x[2],y[2],z[2],0,-V3.dot(x,EYE),-V3.dot(y,EYE),-V3.dot(z,EYE),1];
  const t=1/Math.tan(CURFOV/2),n=0.02,f=60,pr=[t/asp,0,0,0,0,t,0,0,0,0,(f+n)/(n-f),-1,0,0,2*f*n/(n-f),0];
  VP=new Float32Array(16);for(let i=0;i<4;i++)for(let j=0;j<4;j++){let s=0;for(let k=0;k<4;k++)s+=pr[k*4+i]*view[j*4+k];VP[j*4+i]=s;}}
function project(p){const v=[0,0,0,0];for(let i=0;i<4;i++)v[i]=VP[i]*p[0]+VP[4+i]*p[1]+VP[8+i]*p[2]+VP[12+i];return [(v[0]/v[3]*0.5+0.5)*cv.clientWidth,(0.5-v[1]/v[3]*0.5)*cv.clientHeight,v[3]];}
// ---------------------------------------------------------------- the measured forest (reused law, from each line's OWN observed spectrum)
// A at full resolution over the window -> 160 display bins (0.25 Mpc/h) by box averaging; ribbon: hairline where clear, wider and
// darker with absorption, local diffusion at strong absorbers; edges not mirror images.
const pForest=prog(`#version 300 es
in vec4 aP;in vec2 aU;uniform mat4 uVP;uniform float uFace;uniform vec4 uRect;out vec2 vU;out float vA;out vec2 vT;
void main(){vU=aU;vA=aP.w;vT=vec2((aP.x-uRect.x)/uRect.z,(uRect.y-aP.y)/uRect.w);if(uFace>0.5)gl_Position=vec4(vT.x*2.0-1.0,1.0-vT.y*2.0,0.0,1.0);else gl_Position=uVP*vec4(aP.xyz,1.0);}`,`#version 300 es
precision highp float;in vec2 vU;in float vA;in vec2 vT;uniform float uMode;uniform sampler2D uF1;out vec4 o;
float h(vec2 p){return fract(sin(dot(p,vec2(12.9898,78.233)))*43758.5453);}
float n2(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);return mix(mix(h(i),h(i+vec2(1,0)),f.x),mix(h(i+vec2(0,1)),h(i+vec2(1,1)),f.x),f.y);}
void main(){float x=abs(vU.x);
 if(uMode>3.5){o=vec4(vec3(0.03,0.03,0.04),vA*smoothstep(0.72,0.86,x)*(1.0-smoothstep(0.93,1.0,x)));return;}   // set aside: the measured stretch stays, drawn lifted (outline)
 if(uMode>2.5){float s=texture(uF1,vT).r;o=vec4(0.94,0.93,0.91,vA*(1.0-smoothstep(0.6,1.0,x))*smoothstep(0.2,1.0,s));return;}   // paper cut: keeps the line legible inside a wall
 float body=uMode>1.5?pow(max(0.0,1.0-x),2.0):uMode>0.5?1.0-smoothstep(0.55,1.0,x):1.0-smoothstep(0.80,1.0,x);
 float fib=uMode>1.5?0.45+0.55*n2(vec2(vU.x*4.0+vU.y*0.3,vU.y*0.35)):uMode>0.5?1.0:0.86+0.14*n2(vec2(vU.x*9.0,vU.y*0.6));
 o=vec4(vec3(0.03,0.03,0.04),vA*body*fib);}`);
// The ribbon is drawn from the line's observed spectrum at FULL resolution (463 px of 10 km/s over the 40 Mpc/h window): its own forest.
const PXU=M.px_Mpc_h*0.1,pxY=p=>STOP-(p+0.5)*PXU;
const FOREST={},sst=(a,b,x)=>{const t=Math.min(1,Math.max(0,(x-a)/(b-a)));return t*t*(3-2*t);};
function broadOf(A,sig){const n=A.length,R=Math.ceil(3*sig),k=[];for(let d=-R;d<=R;d++)k.push(Math.exp(-0.5*d*d/(sig*sig)));const out=new Float32Array(n);
  for(let p=0;p<n;p++){let s_=0,ws=0;for(let d=-R;d<=R;d++){const j=p+d;if(j<0||j>=n)continue;s_+=k[d+R]*A[j];ws+=k[d+R];}out[p]=sst(0.60,0.95,s_/ws);}return out;}
function forestA(r){const k=ROWIDX(r),A=new Float32Array(NPXW);for(let p=0;p<NPXW;p++)A[p]=1-FLUXF[k*NPXF+X0P+p]/255;return A;}
// one law for every line; given only A[0..n) it uses nothing beyond n (so during a ride, n = the pixels the drop has passed)
function forestLaw(A,r){const n=A.length,As=new Float32Array(n),w=new Float32Array(n),a=new Float32Array(n),hb=new Float32Array(n);
    for(let p=0;p<n;p++)As[p]=0.25*A[Math.max(0,p-1)]+0.5*A[p]+0.25*A[Math.min(n-1,p+1)];
    const feats=[];for(let p=1;p<n-1;p++)if(As[p]>=0.55&&As[p]>=As[p-1]&&As[p]>=As[p+1]){const L=feats[feats.length-1];if(L&&p-L.p<=2){if(As[p]>As[L.p])L.p=p;}else feats.push({p});}
    for(const f of feats)f.s=Math.min(1,(As[f.p]-0.4)/0.5);
    for(let p=0;p<n;p++){const v=As[p];let bl=0;for(const f of feats){const d=p-f.p;if(Math.abs(d)<12)bl+=f.s*Math.exp(-d*d/(2*3.0*3.0));}
      w[p]=0.0016+0.028*Math.pow(Math.max(0,v),1.15)+0.003*Math.min(1,bl);a[p]=Math.min(0.97,0.40+0.57*Math.pow(Math.max(0,v),0.8)+0.15*Math.min(1,bl));hb[p]=Math.min(1,bl);}
    // broad absorption: the same law for every line, from its own spectrum smoothed over 1 Mpc/h. Long saturated stretches (a damped
    // absorber's core and wings, or a saturated wall) read wider and keep a halo; narrow forest absorbers stay knots. No cap equalises them.
    const bw=broadOf(A,1.0/M.px_Mpc_h);for(let p=0;p<n;p++){w[p]+=0.020*bw[p];hb[p]=Math.max(hb[p],bw[p]);}
    // edges: a slight, smooth irregularity (hand-inked), not a serrated leaf margin
    const sL=new Float32Array(n),sR=new Float32Array(n),h1=x=>{const y=Math.sin((r*97+x)*12.9898)*43758.5453;return y-Math.floor(y);};
    for(let p=0;p<n;p++){sL[p]=0.5*h1(2*p)+0.25*h1(2*p-2)+0.25*h1(2*p+2);sR[p]=0.5*h1(2*p+1)+0.25*h1(2*p-1)+0.25*h1(2*p+3);sL[p]=1+0.12*(sL[p]-0.5);sR[p]=1+0.12*(sR[p]-0.5);}
    return {A,As,w,a,hb,sL,sR,feats,bw};}
function forestOf(r){if(!FOREST[r])FOREST[r]=forestLaw(forestA(r),r);return FOREST[r];}
const PART_TAIL=40,PART_WIN=120;
function partialForest(r,np){const full=forestOf(r),p0=Math.max(0,np-PART_WIN),tail=forestLaw(forestA(r).subarray(p0,np),r),cut=Math.max(0,np-PART_TAIL),out={};
  for(const key of ['w','a','hb','bw','sL','sR']){const v=new Float32Array(np);v.set(full[key].subarray(0,cut));for(let p=cut;p<np;p++)v[p]=key==='sL'||key==='sR'?full[key][p]:tail[key][p-p0];out[key]=v;}return out;}
// the damped stretch of a set-aside line, in this tile's pixels (sample i of 20 km/s = pixels 2i, 2i+1)
function liftAmount(r){if(L.phase==='revise'&&L.revise&&L.revise.r===r){const k=ease((L.t-L.tp)/TIM.lift);return L.revise.toAside?k:1-k;}return ST.aside.has(r)?1:0;}
function asideRange(r){if(liftAmount(r)<=0)return null;const d=absorberOf(r);if(!d)return null;const p0=Math.max(0,2*d.i0-X0P),p1=Math.min(NPXW-1,2*d.i1+1-X0P);return p1>=p0?[p0,p1]:null;}
const fP=buf(new Float32Array(NPXW*2*4),gl.DYNAMIC_DRAW),fU=buf(new Float32Array(NPXW*2*2),gl.DYNAMIC_DRAW);
function drawForest(rows,opts){opts=opts||{};const face=!!opts.face,pxu=face?opts.pxu:cv.height/(2*Math.tan(CURFOV/2)),y1=opts.upto==null?NPXW:opts.upto;
  gl.useProgram(pForest.p);gl.uniformMatrix4fv(pForest.u.uVP,false,VP);gl.uniform1f(pForest.u.uFace,face?1:0);gl.uniform4f(pForest.u.uRect,SX0,STOP,SW,SH);gl.depthMask(false);
  for(const r of rows){const np=opts.partial&&opts.partial.r===r?opts.partial.p:NPXW,F=np<NPXW?partialForest(r,np):forestOf(r),x0=rowX(r),fade=opts.fade&&opts.fade.r===r?opts.fade.a:1;if(fade<0.01||np<2)continue;
    const ar=np<NPXW?null:asideRange(r),la=ar?liftAmount(r):0,lift=p=>ar&&p>=ar[0]&&p<=ar[1]?la:0;   // never during a ride
    gl.uniform1i(pForest.u.uF1,1);gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,TXF[1]);
    for(const mode of (ar?[2,3,0,1,4]:[2,3,0,1])){gl.uniform1f(pForest.u.uMode,mode);const P=new Float32Array(np*8),U=new Float32Array(np*4);
      for(let p=0;p<np;p++){const y=pxY(p),z0=face?0.004:LINEZ(),dist=face?1:Math.hypot(x0-EYE[0],y-EYE[1],z0-EYE[2]),px1=dist/pxu,lf=lift(p);
        const w=mode===1?0.75*px1:mode===2?0.004+0.022*F.hb[p]+0.03*F.bw[p]:mode===3?Math.max(F.w[p],0.75*px1)+1.6*px1:Math.max(F.w[p],0.75*px1),
          al=(mode===1?Math.min(0.97,0.70+0.3*F.a[p])*(1-0.55*lf):mode===2?(0.10*F.hb[p]+0.10*F.bw[p])*(1-lf):mode===3?0.55*(1-lf):mode===4?0.85*lf:0.92*F.a[p]*(1-0.80*lf))*fade;
        const eL=mode===1||mode===3?1:F.sL[p],eR=mode===1||mode===3?1:F.sR[p];
        P.set([x0-w*eL,y,z0,al,x0+w*eR,y,z0,al],p*8);U.set([-1,p,1,p],p*4);}
      gl.bindBuffer(gl.ARRAY_BUFFER,fP);gl.bufferSubData(gl.ARRAY_BUFFER,0,P);attrib(pForest,'aP',fP,4);gl.bindBuffer(gl.ARRAY_BUFFER,fU);gl.bufferSubData(gl.ARRAY_BUFFER,0,U);attrib(pForest,'aU',fU,2);
      gl.drawArrays(gl.TRIANGLE_STRIP,0,np*2);}}
  gl.depthMask(true);}
// ---------------------------------------------------------------- state and render
const ST={obs:[],cur:null,prev:null,mix:1,front:null,t:0,hover:-1,aside:new Set((Q_.get('aside')||'').split(',').filter(x=>x).map(Number))};
function setObs(rows){ST.obs=rows.slice();ST.cur=posterior(ST.obs,ST.aside);ST.prev=ST.cur;uploadField(TXF[0],ST.cur);uploadField(TXF[1],ST.cur);ST.mix=1;ST.front=null;}
const REPI=()=>REPV==='A'?0:REPV==='B'?1:2,LINEZ=()=>REPV==='A'?0.004:0.06;
const pBg=prog(`#version 300 es
in vec2 aQ;out vec2 vQ;void main(){vQ=aQ;gl_Position=vec4(aQ,0.999,1.0);}`,`#version 300 es
precision highp float;in vec2 vQ;uniform float uRep;out vec4 o;void main(){float r=length(vQ*vec2(0.75,1.0));
 vec3 c=uRep>0.5?vec3(0.20,0.20,0.215)*(1.0-0.25*r):vec3(0.50,0.49,0.47)*(1.0-0.20*r*r);o=vec4(c,1.0);}`);
// a flat world-space quad (shadow behind the paper in A; the light box behind the curtain in B/C; the reading plane in C)
const pQuad=prog(`#version 300 es
in vec2 aQ;uniform mat4 uVP;uniform vec3 uO,uU,uV;out vec2 vQ;void main(){vQ=aQ;gl_Position=uVP*vec4(uO+uU*aQ.x+uV*aQ.y,1.0);}`,`#version 300 es
precision highp float;in vec2 vQ;uniform vec4 uC;uniform float uSoft;out vec4 o;void main(){vec2 d=abs(vQ);float e=max(d.x,d.y);o=vec4(uC.rgb,uC.a*(1.0-smoothstep(1.0-uSoft,1.0,e)));}`);
function drawQuad(O,U,V,rgba,soft){gl.useProgram(pQuad.p);gl.uniformMatrix4fv(pQuad.u.uVP,false,VP);gl.uniform3fv(pQuad.u.uO,O);gl.uniform3fv(pQuad.u.uU,U);gl.uniform3fv(pQuad.u.uV,V);gl.uniform4f(pQuad.u.uC,...rgba);gl.uniform1f(pQuad.u.uSoft,soft||0.001);
  attrib(pQuad,'aQ',qBuf,2);gl.drawArrays(gl.TRIANGLE_STRIP,0,4);}
const HOOK={post:null,cutY:null};
function drawSheet(face,swatch,nofield){gl.useProgram(pSheet.p);gl.uniformMatrix4fv(pSheet.u.uVP,false,VP);gl.uniform1f(pSheet.u.uFace,face?1:0);gl.uniform4f(pSheet.u.uRect,SX0,STOP,SW,SH);
  gl.uniform1i(pSheet.u.uF0,0);gl.uniform1i(pSheet.u.uF1,1);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,TXF[0]);gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,TXF[1]);
  const fr=ST.front;gl.uniform1f(pSheet.u.uMix,ST.mix);gl.uniform1f(pSheet.u.uFront,fr?1:0);gl.uniform1f(pSheet.u.uFrontR,fr?fr.R:0);gl.uniform1f(pSheet.u.uFrontU,fr?fr.u:0);gl.uniform1f(pSheet.u.uFrontW,fr?fr.W:1);
  gl.uniform1f(pSheet.u.uSwatch,swatch?1:0);gl.uniform1f(pSheet.u.uNoField,nofield?1:0);gl.uniform1f(pSheet.u.uLayer,face?0:(HOOK.layer||0));gl.uniform1f(pSheet.u.uRep,face?(HOOK.faceMat||0):REPI());gl.uniform1f(pSheet.u.uAlpha,1);gl.uniform1f(pSheet.u.uTime,ST.t);
  gl.uniform1f(pSheet.u.uCut,!face&&REPV==='C'&&HOOK.cutY!=null?1:0);gl.uniform1f(pSheet.u.uCutY,HOOK.cutY==null?0:HOOK.cutY);gl.uniform2f(pSheet.u.uRes,cv.width,cv.height);
  attrib(pSheet,'aQ',qBuf,2);gl.drawArrays(gl.TRIANGLE_STRIP,0,4);}
// ---------------------------------------------------------------- the ride, on this sheet's 160 LOS bins of 0.25 Mpc/h
// Causal (nothing ahead of the drop is known): smoothing looks back only; the pace is fixed per bin (slower through absorption), never
// normalised by the whole line (which would reveal its total absorption); an absorber is confirmed, and sounds, when the drop leaves it;
// a clear stretch is confirmed once it has lasted 5 bins; the broad term trails the drop.
const RIDE_BASE=0.0345;   // s per 0.25 Mpc/h bin at A = 0: a typical 40 Mpc/h ride takes about 5-6 s
function rideOf(r){const k=ROWIDX(r),W=WINB,Aobs=new Float32Array(W),As=new Float32Array(W),dt=new Float32Array(W),tb=new Float32Array(W);
  for(let b=0;b<W;b++)Aobs[b]=Math.min(1,Math.max(0,OBS[k*N+X0B+b]*M.obs_scale+M.Abar));
  for(let b=0;b<W;b++)As[b]=b?0.25*Aobs[b-1]+0.75*Aobs[b]:Aobs[b];
  let acc=0;for(let b=0;b<W;b++){dt[b]=RIDE_BASE*(0.7+0.8*As[b]);acc+=dt[b];tb[b]=acc;}
  const feats=[];for(let b=1;b<W;b++){const nx=b<W-1?As[b+1]:-1;if(As[b]>=0.55&&As[b]>=As[b-1]&&As[b]>=nx){if(feats.length&&b-feats[feats.length-1].b<=3)continue;
      feats.push({b,A:As[b],s:Math.min(1,(As[b]-0.4)/0.5),t:b<W-1?tb[b]:tb[b]-dt[b]/2,depth:(b+0.5)*0.25});}}
  const clears=[];for(let b=0;b<W;){if(As[b]<0.12){let e=b;while(e<W&&As[e]<0.12)e++;if(e-b>=5)clears.push({b0:b,b1:e-1,t:tb[b+4]-dt[b+4],depth:(b+0.5)*0.25});b=e;}else b++;}
  const bw=new Float32Array(W);{const sg=4,R_=12;for(let b=0;b<W;b++){let s_=0,ws=0;for(let d=0;d<=R_&&b-d>=0;d++){const kk=Math.exp(-0.5*d*d/(sg*sg));s_+=kk*As[b-d];ws+=kk;}bw[b]=sst(0.60,0.95,s_/ws);}}
  return {r,k,Aobs,As,dt,tb,T:acc,feats,clears,bw};}
function uAt(R,tau){if(tau<=0)return 0;for(let b=0;b<WINB;b++)if(tau<R.tb[b])return (b+1-(R.tb[b]-tau)/R.dt[b])/WINB;return 1;}
// ---- the drop and its wake (reused shaders; the wake lies in the sheet along the line and is drawn only for bins the drop has passed)
const pDrop=prog(`#version 300 es
in vec4 aP;uniform mat4 uVP;uniform float uScale;void main(){gl_Position=uVP*vec4(aP.xyz,1.0);gl_PointSize=uScale*aP.w/gl_Position.w;}`,`#version 300 es
precision highp float;uniform float uLift,uDark;out vec4 o;void main(){vec2 d=gl_PointCoord*2.0-1.0;
 vec2 q=d-vec2(0.0,0.38);float head=length(q)/0.56,w=mix(0.03,0.50,clamp((d.y+0.95)/1.3,0.0,1.0)),tail=d.y<0.38&&d.y>-0.95?abs(d.x)/w:9.0,f=min(head,tail);if(f>1.0)discard;
 float hl=exp(-pow(length(d-vec2(-0.16,0.24))/0.13,2.0));vec3 c=vec3(0.06,0.06,0.07)*(1.0-0.2*uDark);c=mix(c,vec3(0.95),hl*0.55+uLift*0.25);o=vec4(c,1.0-smoothstep(0.85,1.0,f));}`);
const pWake=prog(`#version 300 es
in vec4 aP;in vec2 aU;uniform mat4 uVP;out vec2 vU;out float vA;void main(){gl_Position=uVP*vec4(aP.xyz,1.0);vU=aU;vA=aP.w;}`,`#version 300 es
precision highp float;in vec2 vU;in float vA;out vec4 o;
float h(vec2 p){return fract(sin(dot(p,vec2(12.9898,78.233)))*43758.5453);}
float n2(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);return mix(mix(h(i),h(i+vec2(1,0)),f.x),mix(h(i+vec2(0,1)),h(i+vec2(1,1)),f.x),f.y);}
void main(){float feather=1.0-pow(abs(vU.x),1.5);
 float streak=0.40+0.60*n2(vec2(vU.x*7.0,vU.y*0.45))*n2(vec2(vU.x*16.0+3.0,vU.y*0.9));o=vec4(vec3(0.07,0.07,0.08),vA*feather*streak);}`);
const wakeP=buf(new Float32Array(WINB*4*2*4),gl.DYNAMIC_DRAW),wakeU=buf(new Float32Array(WINB*4*2*2),gl.DYNAMIC_DRAW);
// one wake law (widths scaled to this sheet: lines are 0.15 units apart): bin b's width and darkness at time T
function wakeAt(R,T,b){const age=Math.max(0,T-(R.tb[b]-R.dt[b]/2)),A=R.As[b];let w=(0.006+0.040*Math.pow(A,1.5))*(1+0.6*Math.min(1,age/1.5))+0.020*R.bw[b],al=(0.16+0.70*A)*Math.exp(-age/3.5)+0.10*A+0.30*R.bw[b];
  for(const f of R.feats){const af=T-f.t;if(af<0)continue;const x=b-f.b,sp=Math.min(1,af),c=Math.exp(-Math.pow(x+2.0+2.4*sp,2)/(2*Math.pow(2.4+1.6*sp,2)));
    w+=f.s*(0.03+0.065*(1-Math.exp(-af/0.4)))*c;al+=f.s*0.60*c*Math.exp(-af/1.8);}return [w,Math.min(0.9,al)];}
function drawWake(R,T,u,outF){if(!R||outF<0.01)return 0;const x0=rowX(R.r),z=LINEZ()+0.002,cb=Math.min(WINB-1,Math.floor(u*WINB)),Wb=[],Ab=[],Yb=[];let maxBin=-1;
  for(let b=0;b<=cb;b++){const yb=b<cb?depthY(b):STOP-SH*u;maxBin=b;const [w,a]=wakeAt(R,T,b);Wb.push(w);Ab.push(a*outF);Yb.push(yb);}
  const sm=(X,i)=>0.25*X[Math.max(0,i-1)]+0.5*X[i]+0.25*X[Math.min(X.length-1,i+1)],Ws=Wb.map((_,i)=>sm(Wb,i)),As_=Ab.map((_,i)=>sm(Ab,i)),P=[],U=[];
  for(let i=0;i<Ws.length;i++){const n=i<Ws.length-1?2:1;for(let j=0;j<n;j++){const f=j/n,i2=Math.min(Ws.length-1,i+1),w=Ws[i]+(Ws[i2]-Ws[i])*f,a=As_[i]+(As_[i2]-As_[i])*f,yb=Yb[i]+(Yb[i2]-Yb[i])*f;
    const fl=0.22*Math.sin(0.65*(i+f)+0.9*T+0.7*R.r),cx_=0.28*w*Math.sin(0.28*(i+f)-1.7*T);
    for(const sg of [-1,1]){P.push(x0+cx_+sg*w*(1+sg*fl),yb,z,a);U.push(sg,(i+f)*0.5);}}}
  if(P.length<8)return maxBin;gl.useProgram(pWake.p);gl.uniformMatrix4fv(pWake.u.uVP,false,VP);
  gl.bindBuffer(gl.ARRAY_BUFFER,wakeP);gl.bufferSubData(gl.ARRAY_BUFFER,0,new Float32Array(P));attrib(pWake,'aP',wakeP,4);gl.bindBuffer(gl.ARRAY_BUFFER,wakeU);gl.bufferSubData(gl.ARRAY_BUFFER,0,new Float32Array(U));attrib(pWake,'aU',wakeU,2);
  gl.drawArrays(gl.TRIANGLE_STRIP,0,P.length/4);return maxBin;}
function drawDrop(R,u,T){const x0=rowX(R.r),y=STOP-SH*u,A=R.As[Math.min(WINB-1,Math.floor(u*WINB))];gl.useProgram(pDrop.p);gl.uniformMatrix4fv(pDrop.u.uVP,false,VP);gl.uniform1f(pDrop.u.uScale,cv.height*0.085);
  gl.uniform1f(pDrop.u.uLift,A<0.12?1:0);gl.uniform1f(pDrop.u.uDark,A);gl.bindBuffer(gl.ARRAY_BUFFER,dynBuf);gl.bufferSubData(gl.ARRAY_BUFFER,0,new Float32Array([x0,y+0.012,LINEZ()+0.004,1.0]));attrib(pDrop,'aP',dynBuf,4);gl.drawArrays(gl.POINTS,0,1);}
const dynBuf=buf(new Float32Array(8192),gl.DYNAMIC_DRAW);
// ---- candidate sources above the top edge (one per candidate line; geometry only: they never show a candidate's data)
const pLant=prog(`#version 300 es
in vec4 aP;uniform mat4 uVP;void main(){gl_Position=uVP*vec4(aP.xyz,1.0);gl_PointSize=abs(aP.w);}`,`#version 300 es
precision highp float;uniform float uSpent,uAlpha,uHot;out vec4 o;void main(){vec2 d=gl_PointCoord*2.0-1.0;float r=length(d);if(r>1.0)discard;
 if(uSpent>0.5){o=vec4(0.08,0.08,0.09,uAlpha*(1.0-smoothstep(0.40,0.62,r)));return;}float core=1.0-smoothstep(0.16,0.30,r),halo=exp(-pow((r-0.52)/0.10,2.0)),glow=exp(-r*r*3.0);
 o=vec4(mix(vec3(0.12,0.12,0.14),vec3(1.0,0.99,0.95),clamp(core+0.7*glow*(1.0-halo)+uHot*0.3,0.0,1.0)),uAlpha*clamp(core+0.8*halo+(0.30+0.4*uHot)*glow,0.0,1.0));}`);
const srcPos=r=>[rowX(r),STOP+0.085,LINEZ()];
function drawSources(){const px=cv.height/720;gl.useProgram(pLant.p);gl.uniformMatrix4fv(pLant.u.uVP,false,VP);gl.bindBuffer(gl.ARRAY_BUFFER,dynBuf);attrib(pLant,'aP',dynBuf,4);
  for(const r of SURV){const sp=ST.obs.includes(r),hot=r===ST.hover&&L.phase==='choose';const p=srcPos(r);gl.uniform1f(pLant.u.uSpent,sp?1:0);gl.uniform1f(pLant.u.uHot,hot?1:0);
    gl.uniform1f(pLant.u.uAlpha,L.phase==='choose'||sp?1:0.45);gl.bufferSubData(gl.ARRAY_BUFFER,0,new Float32Array([...p,(sp?9:hot?24:17)*px]));gl.drawArrays(gl.POINTS,0,1);}}
// hover guide: where the line WOULD run (a plain dotted column; no data)
const pFlat=prog(`#version 300 es
in vec3 aP;uniform mat4 uVP;void main(){gl_Position=uVP*vec4(aP,1.0);}`,`#version 300 es
precision highp float;uniform vec4 uC;out vec4 o;void main(){o=uC;}`);
function drawFlat(arr,mode,rgba){gl.useProgram(pFlat.p);gl.uniformMatrix4fv(pFlat.u.uVP,false,VP);gl.uniform4f(pFlat.u.uC,...rgba);gl.bindBuffer(gl.ARRAY_BUFFER,dynBuf);gl.bufferSubData(gl.ARRAY_BUFFER,0,new Float32Array(arr));attrib(pFlat,'aP',dynBuf,3);gl.drawArrays(mode,0,arr.length/3);}
function drawGuide(r,a){const x=rowX(r),z=LINEZ()+0.001,arr=[];for(let y=STOP;y>STOP-SH;y-=0.06)arr.push(x,y,z,x,Math.max(STOP-SH,y-0.03),z);drawFlat(arr,gl.LINES,[0.15,0.15,0.17,0.55*a]);}
// ---- C: one restrained transverse reading plane: the inferred profile across the sheet at one depth, drawn as relief out of the sheet
function fieldAt(t,b){const f=ST.front,sN=ST.cur.s[t*WINB+b],cN=cAt(ST.cur,t,b);if(!f)return [sN,cN];const d=Math.abs((t+0.5)/NR-f.u)*NR*0.25,m=Math.min(1,Math.max(0,(f.R-(d-f.W))/f.W));
  return [ST.prev.s[t*WINB+b]*(1-m)+sN*m,cAt(ST.prev,t,b)*(1-m)+cN*m];}
function drawCut(){if(REPV!=='C'||HOOK.cutY==null)return;const y=HOOK.cutY,b=Math.min(WINB-1,Math.max(0,Math.round((STOP-y)/SH*WINB-0.5)));
  drawQuad([0,y,0.02],[SW/2+0.06,0,0],[0,0,0.34],[0.96,0.95,0.92,0.16],0.08);
  const thick=[],thin=[];for(let t=0;t<NR-1;t++){const [s0,c0]=fieldAt(t,b),[s1]=fieldAt(t+1,b),seg=[rowX(LO+t),y,0.02+0.16*s0,rowX(LO+t+1),y,0.02+0.16*s1];(c0>=0.5?thick:thin).push(...seg);}
  drawFlat(thin,gl.LINES,[0.12,0.12,0.14,0.35]);for(const dz of [-0.0015,0,0.0015])drawFlat(thick.map((v,i)=>i%3===1?v+dz:v),gl.LINES,[0.08,0.08,0.10,0.85]);}
// ---- gold: returned where the inferred structure changed (from the new posterior only; secondary)
const pGold=prog(`#version 300 es
in vec4 aP;uniform mat4 uVP;void main(){gl_Position=uVP*vec4(aP.xyz,1.0);gl_PointSize=aP.w;}`,`#version 300 es
precision highp float;uniform float uA;out vec4 o;void main(){vec2 d=gl_PointCoord*2.0-1.0;float r=length(d);if(r>1.0)discard;float core=1.0-smoothstep(0.40,0.55,r),ring=exp(-pow((r-0.62)/0.08,2.0)),glow=exp(-r*r*2.5);
 vec3 c=mix(vec3(0.86,0.64,0.16),vec3(1.0,0.92,0.62),(1.0-smoothstep(0.0,0.35,r))*0.6);c=mix(c,vec3(0.30,0.20,0.04),ring*0.7);o=vec4(c,uA*clamp(core+0.9*ring+0.25*glow,0.0,1.0));}`);
function goldSeeds(prev,cur,r){const t0=r-LO,S=[];   // local maxima of |change| x new support, beyond the line's own cells
  const g=(t,b)=>Math.abs(t-t0)<3?-1:Math.abs(cur.s[t*WINB+b]-prev.s[t*WINB+b])*cAt(cur,t,b);
  for(let t=2;t<NR-2;t+=2)for(let b=3;b<WINB-3;b+=2){if(Math.abs(t-t0)<3)continue;const v=g(t,b);if(v<0.5)continue;let mx=true;for(let dt=-4;dt<=4&&mx;dt+=2)for(let db=-6;db<=6;db+=2)if((dt||db)&&t+dt>=0&&t+dt<NR&&b+db>=0&&b+db<WINB&&g(t+dt,b+db)>v){mx=false;break;}
    if(mx)S.push({t,b,v,p:[rowX(LO+t),depthY(b),LINEZ()+0.003],d:Math.abs(t-t0)*0.25});}
  S.sort((a,b)=>b.v-a.v);return S.slice(0,7).sort((a,b)=>a.d-b.d);}
// ---------------------------------------------------------------- camera poses (per representation)
const ease=x=>{x=Math.min(1,Math.max(0,x));return x*x*(3-2*x);},lerp=(a,b,t)=>a+(b-a)*t,lerp3=(a,b,t)=>[lerp(a[0],b[0],t),lerp(a[1],b[1],t),lerp(a[2],b[2],t)];
function poseLerp(A,B,t){return {eye:lerp3(A.eye,B.eye,t),tgt:lerp3(A.tgt,B.tgt,t),fov:lerp(A.fov,B.fov,t)};}
function overPose(){const k=SW/3,asp=Math.max(0.3,(cv.clientWidth||16)/(cv.clientHeight||9)),zf=(SW/2+0.22)/(Math.tan(0.43)*asp)+0.4;
  return REPV==='A'?{eye:[0.46*k,0.28,Math.max(4.0*k+0.3,zf)],tgt:[0,asp<1?0.20:0.10,0],fov:0.86}:{eye:[2.05*k,0.50,Math.max(3.55*k+0.4,zf)],tgt:[0.10,0.10,0],fov:0.86};}
function emergePose(r){const O=overPose(),x=rowX(r)*0.4,k=0.86;return {eye:[x+(O.eye[0]-O.tgt[0])*k,O.eye[1]*0.9,O.eye[2]*k],tgt:[x,O.tgt[1]-0.05,0],fov:O.fov};}
function divePose(r,u){const x=rowX(r),y=Math.max(STOP-SH*u,STOP-SH+0.42),side=REPV==='A'?0.10:0.34;return {eye:[x+side,y-0.05,LINEZ()+1.20],tgt:[x+side*0.3,y-0.20,LINEZ()],fov:0.80};}
// ---------------------------------------------------------------- audio (reused synthesis: one event list, rendered live or offline for recordings)
const LIVE=MODE==='play';
const AU={ctx:null,out:null,ev:[],clock:0};
function auChain(ctx){const m=ctx.createGain();m.gain.value=1.6;const comp=ctx.createDynamicsCompressor();comp.threshold.value=-14;m.connect(comp);comp.connect(ctx.destination);
  const verb=ctx.createConvolver(),n=Math.floor(ctx.sampleRate*2.6),bf=ctx.createBuffer(2,n,ctx.sampleRate);for(let c=0;c<2;c++){const d=bf.getChannelData(c);let s=12345+c;for(let i=0;i<n;i++){s=(s*16807)%2147483647;d[i]=(s/2147483647*2-1)*Math.pow(1-i/n,2.6);}}
  verb.buffer=bf;const vg=ctx.createGain();vg.gain.value=0.3;verb.connect(vg);vg.connect(m);return {m,verb};}
function nb(ctx,sec){const b=ctx.createBuffer(1,Math.max(1,Math.floor(ctx.sampleRate*sec)),ctx.sampleRate),d=b.getChannelData(0);let s=777;for(let i=0;i<d.length;i++){s=(s*16807)%2147483647;d[i]=s/2147483647*2-1;}return b;}
function send(ctx,o,node,wet){const g=ctx.createGain();node.connect(g);g.connect(o.m);const w=ctx.createGain();w.gain.value=wet;node.connect(w);w.connect(o.verb);}
const NOTE=k=>196*Math.pow(2,[0,2,4,7,9,12,14,16,19,21][k%10]/12+Math.floor(k/10));
const SYN={
  launch(ctx,o,t){const s=ctx.createOscillator(),g=ctx.createGain();s.type='sine';s.frequency.setValueAtTime(880,t);s.frequency.exponentialRampToValueAtTime(330,t+0.18);
    g.gain.setValueAtTime(0.0001,t);g.gain.exponentialRampToValueAtTime(0.10,t+0.01);g.gain.exponentialRampToValueAtTime(0.0001,t+0.35);s.connect(g);send(ctx,o,g,0.4);s.start(t);s.stop(t+0.4);},
  dive(ctx,o,t,a){const T=a.T,n=ctx.createBufferSource();n.buffer=nb(ctx,T+2);const bp=ctx.createBiquadFilter();bp.type='bandpass';bp.Q.value=0.7;const g=ctx.createGain();
    const v=ctx.createBufferSource();v.buffer=nb(ctx,T+2);const vbp=ctx.createBiquadFilter();vbp.type='bandpass';vbp.Q.value=2.2;const vg=ctx.createGain();
    const lo=ctx.createOscillator();lo.type='triangle';const lg=ctx.createGain();lo.connect(lg);
    g.gain.setValueAtTime(0.0001,t);vg.gain.setValueAtTime(0.0001,t);lg.gain.setValueAtTime(0.0001,t);let prev=0;
    for(let b=0;b<a.As.length;b++){const tt=t+prev,A=a.As[b],spd=1/a.dt[b];prev=a.tb[b];
      bp.frequency.setValueAtTime(400+1200*Math.min(1,spd*0.03),tt);g.gain.linearRampToValueAtTime(0.03+0.03*Math.min(1,spd*0.03),tt);
      vbp.frequency.setValueAtTime(420+1600*A,tt);vg.gain.linearRampToValueAtTime(0.0006+0.06*A*A,tt);lo.frequency.setValueAtTime(NOTE(0)/2*(1-0.25*A),tt);lg.gain.linearRampToValueAtTime(0.002+0.05*A,tt);}
    for(const gg of [g,vg,lg])gg.gain.setTargetAtTime(0.0001,t+T,0.15);
    n.connect(bp);bp.connect(g);send(ctx,o,g,0.3);v.connect(vbp);vbp.connect(vg);send(ctx,o,vg,0.25);send(ctx,o,lg,0.5);n.start(t);v.start(t);lo.start(t);n.stop(t+T+0.7);v.stop(t+T+0.7);lo.stop(t+T+0.7);},
  hit(ctx,o,t,a){const c=ctx.createOscillator(),m=ctx.createOscillator(),mg=ctx.createGain(),g=ctx.createGain(),f=NOTE(2)/2*(1+0.3*(1-a.s));c.frequency.value=f;m.frequency.value=f*1.41;
    mg.gain.setValueAtTime(f*2.4,t);mg.gain.exponentialRampToValueAtTime(1,t+1.2);m.connect(mg);mg.connect(c.frequency);c.connect(g);g.gain.setValueAtTime(0.0001,t);g.gain.exponentialRampToValueAtTime(0.04+0.10*a.s,t+0.006);
    g.gain.exponentialRampToValueAtTime(0.0001,t+1.4);send(ctx,o,g,0.5);c.start(t);m.start(t);c.stop(t+1.6);m.stop(t+1.6);
    const th=ctx.createOscillator(),tg=ctx.createGain();th.frequency.setValueAtTime(90,t);th.frequency.exponentialRampToValueAtTime(48,t+0.25);tg.gain.setValueAtTime(0.0001,t);tg.gain.exponentialRampToValueAtTime(0.06*a.s+0.001,t+0.005);tg.gain.exponentialRampToValueAtTime(0.0001,t+0.3);th.connect(tg);send(ctx,o,tg,0.1);th.start(t);th.stop(t+0.35);},
  clear(ctx,o,t){const s=ctx.createOscillator(),g=ctx.createGain();s.frequency.value=NOTE(7);g.gain.setValueAtTime(0.0001,t);g.gain.exponentialRampToValueAtTime(0.018,t+0.05);g.gain.exponentialRampToValueAtTime(0.0001,t+0.9);s.connect(g);send(ctx,o,g,0.8);s.start(t);s.stop(t+1.0);},
  emerge(ctx,o,t){const n=ctx.createBufferSource();n.buffer=nb(ctx,1.4);const bp=ctx.createBiquadFilter();bp.type='bandpass';bp.Q.value=0.8;bp.frequency.setValueAtTime(300,t);bp.frequency.exponentialRampToValueAtTime(2400,t+1.0);
    const g=ctx.createGain();g.gain.setValueAtTime(0.0001,t);g.gain.exponentialRampToValueAtTime(0.05,t+0.35);g.gain.exponentialRampToValueAtTime(0.0001,t+1.2);n.connect(bp);bp.connect(g);send(ctx,o,g,0.4);n.start(t);n.stop(t+1.4);},
  brush(ctx,o,t,a){const d=a.d||0.62,n=ctx.createBufferSource();n.buffer=nb(ctx,d+0.3);const hp=ctx.createBiquadFilter();hp.type='highpass';hp.frequency.value=1400;const bp=ctx.createBiquadFilter();bp.type='bandpass';bp.Q.value=1.3;
    bp.frequency.setValueAtTime(2200,t);bp.frequency.linearRampToValueAtTime(3800,t+d);const am=ctx.createGain(),g=ctx.createGain();am.gain.setValueAtTime(0.0,t);
    let sd=31;for(let i=0,tt=t;tt<t+d;i++){sd=(sd*16807)%2147483647;const r=sd/2147483647;tt+=0.018+0.03*r;am.gain.setValueAtTime(0.25+0.75*r,tt);}
    g.gain.setValueAtTime(0.0001,t);g.gain.exponentialRampToValueAtTime(0.05,t+0.16);g.gain.setValueAtTime(0.05,t+d*0.6);g.gain.exponentialRampToValueAtTime(0.0001,t+d);
    n.connect(hp);hp.connect(bp);bp.connect(am);am.connect(g);send(ctx,o,g,0.18);n.start(t);n.stop(t+d+0.2);},
  gold(ctx,o,t,a){const k=a.amp||1;for(const [p,amp,d] of [[1,0.10,1.4],[2.76,0.035,0.8],[5.4,0.015,0.4]]){const s=ctx.createOscillator(),g=ctx.createGain();s.frequency.value=NOTE(9+a.shot)*p;g.gain.setValueAtTime(0.0001,t);g.gain.exponentialRampToValueAtTime(amp*k,t+0.003);g.gain.exponentialRampToValueAtTime(0.0001,t+d);s.connect(g);send(ctx,o,g,0.55);s.start(t);s.stop(t+d+0.1);}},
};
function au(type,t,a){AU.ev.push({type,t,a});if(LIVE&&AU.ctx){try{SYN[type](AU.ctx,AU.out,AU.ctx.currentTime+Math.max(0,t-AU.clock),a||{});}catch(e){}}}
function audioInit(){if(AU.ctx&&AU.ctx.state==='suspended'){try{AU.ctx.resume();}catch(e){}}if(AU.ctx||!LIVE)return;const C=window.AudioContext||window.webkitAudioContext;if(!C)return;AU.ctx=new C();AU.out=auChain(AU.ctx);}
async function audioWav(dur){const sr=44100,ctx=new OfflineAudioContext(2,Math.ceil(sr*(dur+2.5)),sr),o=auChain(ctx);for(const e of AU.ev)SYN[e.type](ctx,o,e.t,e.a||{});
  const bf=await ctx.startRendering(),n=bf.length,ch=[bf.getChannelData(0),bf.getChannelData(1)],dv=new DataView(new ArrayBuffer(44+n*4));
  const w=(o,s)=>{for(let i=0;i<s.length;i++)dv.setUint8(o+i,s.charCodeAt(i));};w(0,'RIFF');dv.setUint32(4,36+n*4,true);w(8,'WAVE');w(12,'fmt ');dv.setUint32(16,16,true);dv.setUint16(20,1,true);dv.setUint16(22,2,true);
  dv.setUint32(24,sr,true);dv.setUint32(28,sr*4,true);dv.setUint16(32,4,true);dv.setUint16(34,16,true);w(36,'data');dv.setUint32(40,n*4,true);
  for(let i=0;i<n;i++)for(let c=0;c<2;c++)dv.setInt16(44+(i*2+c)*2,Math.max(-1,Math.min(1,ch[c][i]))*32767,true);
  const u=new Uint8Array(dv.buffer);let s='';for(let i=0;i<u.length;i+=32768)s+=String.fromCharCode.apply(null,u.subarray(i,i+32768));return btoa(s);}
// ---------------------------------------------------------------- the loop: CHOOSE -> DIVE (the measured line) -> EMERGE (the field revises) -> CHOOSE ... -> DONE (hold) -> DEBRIEF
// One player action besides choosing: press-and-hold a landed ribbon to set aside its damped stretch, or to restore it.
const BUDGET=Math.max(1,Math.min(20,+(Q_.get('budget')||10)));   // 25 candidate sightlines, 10 drops
const TIM={launch:0.8,dry:0.6,pull:1.3,front:2.8,gold:1.0,settle:1.0,hold:0.9,lift:0.6,endHold:7.0};
const L={phase:'choose',t:0,tp:0,shot:0,R:null,r:-1,log:[],gold:0,seeds:[],fromPose:null,dvMax:0,hold:null,revise:null,nPost:0,msgUntil:0};
function log(ev,o){L.log.push(Object.assign({ev,t:+L.t.toFixed(3),shot:L.shot},o||{}));}
const HINTS={};function hint(k,txt,secs){if(HINTS[k]||MODE==='still')return;HINTS[k]=1;say(txt);L.msgUntil=L.t+(secs||6);}
const TOUCH=matchMedia('(pointer:coarse)').matches||Q_.get('touch')==='1';
function start(r){if(L.phase!=='choose'||ST.obs.includes(r)||!SURV.includes(r)||L.shot>=BUDGET)return false;
  L.r=r;L.R=rideOf(r);L.phase='dive';L.tp=L.t;L.fromPose=JSON.parse(JSON.stringify(CAM.pose));ST.hover=-1;L.hold=null;
  au('launch',L.t);au('dive',L.t+TIM.launch,{T:L.R.T,As:Array.from(L.R.As),dt:Array.from(L.R.dt),tb:Array.from(L.R.tb)});
  for(const f of L.R.feats)au('hit',L.t+TIM.launch+f.t,{s:f.s});for(const c of L.R.clears)au('clear',L.t+TIM.launch+c.t);
  log('choose',{row:r,T:+L.R.T.toFixed(2)});say('');hint('ride','The ride is this sightline’s measured Lyα forest: slower and darker where the gas absorbs.',5.5);return true;}
function beginEmerge(){const prev=ST.cur,t0=performance.now(),nxt=posterior(ST.obs.concat([L.r]),ST.aside,MODE!=='still'),ms=performance.now()-t0;L.cjob=nxt.job||null;L.cjobT0=performance.now();L.nPost++;   // the posterior is updated only now: the observation is complete
  ST.prev=prev;ST.cur=nxt;ST.obs=ST.obs.concat([L.r]);uploadField(TXF[0],prev);if(!L.cjob)uploadField(TXF[1],nxt);
  const u=(L.r-LO+0.5)/NR,Rmax=Math.max(u,1-u)*NR*0.25+4.0;ST.front={u,R:0,W:3.5,Rmax};L.seeds=L.cjob?[]:goldSeeds(prev,nxt,L.r);
  let dv=0;if(!L.cjob)for(let i=0;i<nxt.s.length;i++)dv=Math.max(dv,Math.abs(nxt.s[i]-prev.s[i]));L.dvMax=dv;
  L.phase='emerge';L.tp=L.t;L.fromPose=JSON.parse(JSON.stringify(CAM.pose));HOOK.forestFade={r:L.r,a:0.8};
  au('brush',L.t,{d:TIM.dry});au('emerge',L.t+TIM.dry+TIM.pull*0.6);
  if(REPV==='C'){let best=-1,bb=WINB>>1;for(let b=0;b<WINB;b++){let v=0;for(let t=0;t<NR;t++)v+=Math.abs(nxt.s[t*WINB+b]-prev.s[t*WINB+b])*cAt(nxt,t,b);if(v>best){best=v;bb=b;}}L.cutTo=depthY(bb);}
  log('emerge',{row:L.r,solve_ms:+ms.toFixed(1),dv_max:+dv.toFixed(3),gold:L.seeds.length,n:ST.obs.length,masked:nxt.masked||[]});}
function endEmerge(){uploadField(TXF[0],ST.cur);ST.front=null;ST.prev=ST.cur;HOOK.forestFade=null;L.gold+=L.seeds.length;L.shot++;L.R=null;
  L.phase=L.shot>=BUDGET?'done':'choose';L.tp=L.t;L.fromPose=JSON.parse(JSON.stringify(CAM.pose));log(L.phase==='done'?'done':'ready',{n:ST.obs.length,gold:L.gold});budgetUI();
  if(L.shot===1)hint('sheet','The sheet shows the large-scale structure inferred between your sightlines. Bare paper is still unknown.',6.5);
  if(L.shot===2)hint('hold','Press and hold a landed line to set it aside. Hold it again to restore it.',6.5);
  if(L.phase==='done')say('');}
// ---- set aside / restore: the damped-absorber rule on the line's own observed spectrum; exact masked conditioning; no gold, no verdict
function holdComplete(r){if(!ST.obs.includes(r))return false;const d=absorberOf(r);
  if(!d){log('hold_noop',{row:r});return false;}   // no damped stretch on this line: nothing happens
  const toAside=!ST.aside.has(r);if(toAside)ST.aside.add(r);else ST.aside.delete(r);
  log(toAside?'aside':'restore',{row:r,mask_samples:[d.i0,d.i1],mask_Mpc:[+d.Mpc[0].toFixed(2),+d.Mpc[1].toFixed(2)],n_aside:ST.aside.size});
  beginRevise(r,toAside);return true;}
function beginRevise(r,toAside){const from=L.phase==='debrief'||L.phase==='done'?'done':'choose';if(L.phase==='debrief')hideDebrief();
  const prev=ST.cur,t0=performance.now(),nxt=posterior(ST.obs,ST.aside,MODE!=='still'),ms=performance.now()-t0;L.nPost++;
  L.cjob=nxt.job||null;L.cjobT0=performance.now();ST.prev=prev;ST.cur=nxt;uploadField(TXF[0],prev);if(!L.cjob)uploadField(TXF[1],nxt);
  const u=(r-LO+0.5)/NR,Rmax=Math.max(u,1-u)*NR*0.25+4.0;ST.front={u,R:0,W:3.5,Rmax};
  let dv=0;if(!L.cjob)for(let i=0;i<nxt.s.length;i++)dv=Math.max(dv,Math.abs(nxt.s[i]-prev.s[i]));
  L.revise={r,toAside,from,t0,frontAt:null};L.phase='revise';L.tp=L.t;L.fromPose=JSON.parse(JSON.stringify(CAM.pose));L.R=null;say('');
  au('brush',L.t,{d:TIM.lift});log('revise_begin',{row:r,aside:toAside,mean_ms:+ms.toFixed(1),dv_max:+dv.toFixed(3),masked:nxt.masked||[]});}
function endRevise(){uploadField(TXF[0],ST.cur);ST.front=null;ST.prev=ST.cur;const from=L.revise.from;log('revise_end',{row:L.revise.r});L.revise=null;
  L.phase=from;L.tp=L.t;L.fromPose=JSON.parse(JSON.stringify(CAM.pose));}
const ET=()=>TIM.dry+TIM.pull+TIM.front+TIM.gold;
function step(dt){L.t+=dt;ST.t=L.t;const tau=L.t-L.tp;if(L.msgUntil&&L.t>L.msgUntil){L.msgUntil=0;say('');}
  if(L.phase==='choose'||L.phase==='done'||L.phase==='debrief'){const P=L.phase==='choose'?overPose():L.phase==='debrief'?debriefPose():donePose(),k=L.phase==='choose'?TIM.settle:3.0;CAM.pose=L.fromPose?poseLerp(L.fromPose,P,ease(tau/k)):P;if(REPV==='C'&&HOOK.cutY==null)HOOK.cutY=0.4;
    if(L.phase==='done'&&tau>=TIM.endHold&&!L.hold)showDebrief();}
  else if(L.phase==='dive'){const R=L.R,td=tau-TIM.launch,u=uAt(R,td);
    CAM.pose=td<0?poseLerp(L.fromPose,divePose(R.r,0),ease(tau/TIM.launch)):divePose(R.r,u);if(td>=R.T)beginEmerge();}
  else if(L.phase==='emerge'){const R=L.R;
    if(L.cjob){const done=MODE==='play'?L.cjob.run(1e9,7):L.cjob.run(8);if(done){const J=L.cjob;L.cjob=null;uploadField(TXF[1],ST.cur);L.seeds=goldSeeds(ST.prev,ST.cur,L.r);log('support',{row:L.r,ms:+J.ms.toFixed(1),mean_ms:+(J.meanMs||0).toFixed(1),units:J.units,wall_ms:+(performance.now()-L.cjobT0).toFixed(1)});}}
    if(tau<TIM.dry){HOOK.forestFade={r:L.r,a:0.8+0.2*ease(tau/TIM.dry)};CAM.pose=divePose(R.r,1);}
    else{HOOK.forestFade=null;const tp=tau-TIM.dry;CAM.pose=poseLerp(divePose(R.r,1),emergePose(R.r),ease(tp/TIM.pull));
      const tf=tau-TIM.dry-TIM.pull*0.55;if(ST.front&&!L.cjob)ST.front.R=Math.max(0,tf)/TIM.front*ST.front.Rmax;
      if(REPV==='C'&&L.cutTo!=null)HOOK.cutY=lerp(HOOK.cutY,L.cutTo,1-Math.exp(-dt*2.5));
      for(const g of L.seeds){const ta=TIM.dry+TIM.pull*0.55+TIM.front*Math.min(1,(g.d+3.5)/ST.front.Rmax);if(!g.au&&tau>=ta){g.au=1;au('gold',L.t,{shot:L.shot,amp:0.5});}}}
    if(tau>=ET())endEmerge();}
  else if(L.phase==='revise'){const V=L.revise;
    if(L.cjob){const done=MODE==='play'?L.cjob.run(1e9,7):L.cjob.run(8);if(done){const J=L.cjob;L.cjob=null;uploadField(TXF[1],ST.cur);let dv=0;for(let i=0;i<ST.cur.s.length;i++)dv=Math.max(dv,Math.abs(ST.cur.s[i]-ST.prev.s[i]));log('support',{row:V.r,ms:+J.ms.toFixed(1),mean_ms:+(J.meanMs||0).toFixed(1),units:J.units,wall_ms:+(performance.now()-L.cjobT0).toFixed(1),revise:true,dv_max:+dv.toFixed(3)});}}
    CAM.pose=poseLerp(L.fromPose,emergePose(V.r),ease(tau/1.0));
    if(!L.cjob&&V.frontAt==null&&tau>=TIM.lift){V.frontAt=L.t;au('emerge',L.t);log('revise_front',{row:V.r,latency_s:+(L.t-L.tp).toFixed(3),wall_ms:+(performance.now()-V.t0).toFixed(1)});}
    if(V.frontAt!=null){const tf=L.t-V.frontAt;ST.front.R=tf/TIM.front*ST.front.Rmax;if(tf>=TIM.front+0.6)endRevise();}}}
function donePose(){const O=overPose();return {eye:[O.eye[0]*1.18,O.eye[1]+0.15,O.eye[2]*1.04],tgt:O.tgt,fov:O.fov};}
// while the debrief card is open the sheet slides clear of it (left card on wide screens, bottom sheet on narrow ones)
function debriefPose(){const P=donePose(),open=document.getElementById('deb').style.display==='block';if(!open)return P;const narrow=innerWidth<=640,d=narrow?[0,-1.15,0]:[-0.95*SW/3.75,0,0];
  return {eye:V3.add(P.eye,d),tgt:V3.add(P.tgt,d),fov:P.fov};}
function drawGold(){if(L.phase!=='emerge'||!ST.front)return;const tau=L.t-L.tp,px=cv.height/720,home=[SX0+0.05,STOP+0.32,LINEZ()],pts=[];
  for(const g of L.seeds){const ta=TIM.dry+TIM.pull*0.55+TIM.front*Math.min(1,(g.d+3.5)/ST.front.Rmax),tf=ET()-TIM.gold*0.95;if(tau<ta)continue;const al=Math.min(1,(tau-ta)/0.25);
    if(tau<tf){pts.push([g.p,10*px*(0.9+0.15*Math.sin((tau-ta)*9+g.t)),al]);continue;}const a=ease((tau-tf)/(TIM.gold*0.9));if(a>=1)continue;
    const mid=[(g.p[0]+home[0])/2,Math.max(g.p[1],home[1])+0.4,0.2];pts.push([lerp3(lerp3(g.p,mid,a),lerp3(mid,home,a),a),10*px,1]);}
  gl.useProgram(pGold.p);gl.uniformMatrix4fv(pGold.u.uVP,false,VP);gl.bindBuffer(gl.ARRAY_BUFFER,dynBuf);attrib(pGold,'aP',dynBuf,4);
  for(const [q,sz,aa] of pts){gl.uniform1f(pGold.u.uA,0.95*aa);gl.bufferSubData(gl.ARRAY_BUFFER,0,new Float32Array([q[0],q[1],q[2],sz]));gl.drawArrays(gl.POINTS,0,1);}}
function render(){const W_=Math.round(cv.clientWidth*devicePixelRatio),H_=Math.round(cv.clientHeight*devicePixelRatio);if(cv.width!==W_||cv.height!==H_){cv.width=W_;cv.height=H_;}
  gl.viewport(0,0,cv.width,cv.height);camUpdate();gl.clearColor(0.3,0.3,0.3,1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);
  gl.disable(gl.DEPTH_TEST);gl.useProgram(pBg.p);gl.uniform1f(pBg.u.uRep,REPI());attrib(pBg,'aQ',qBuf,2);gl.drawArrays(gl.TRIANGLE_STRIP,0,4);
  gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);
  if(REPV==='A')drawQuad([0.05,-0.07,-0.05],[SW/2+0.05,0,0],[0,SH/2+0.05,0],[0.05,0.05,0.05,0.45],0.06);
  else{for(const [k,a] of [[2.3,0.10],[1.8,0.18]])drawQuad([0,0,-0.55],[SW/2*k,0,0],[0,SH/2*k,0],[1.0,0.975,0.92,a],0.55);drawQuad([0,0,-0.5],[SW/2*1.6,0,0],[0,SH/2*1.45,0],[1.0,0.99,0.96,1.0],0.30);}   // the light behind the curtain is uniform over the whole curtain from any play angle
  drawSheet(false,false);if(REPV==='C')drawCut();
  if(L.phase==='choose'&&ST.hover>=0)drawGuide(ST.hover,1);
  const ff=HOOK.forestFade,lf=HOOK.layer===2?{fade:{a:0.35,all:1}}:{};drawForest(ST.obs.concat(ff&&!ST.obs.includes(ff.r)?[ff.r]:[]),ff?{fade:ff}:lf);
  if(L.R&&L.phase==='dive'){const td=L.t-L.tp-TIM.launch;if(td>0){const u=uAt(L.R,td);drawForest([L.R.r],{partial:{r:L.R.r,p:Math.floor(u*NPXW)},fade:{r:L.R.r,a:0.8}});}}
  if(L.R&&(L.phase==='dive'||L.phase==='emerge')){const tau=L.t-L.tp,td=L.phase==='dive'?tau-TIM.launch:L.R.T,u=L.phase==='dive'?uAt(L.R,td):1,
      outF=L.phase==='dive'?1:1-ease(tau/TIM.dry);L.maxBin=drawWake(L.R,Math.max(0,td),u,outF);if(L.phase==='dive'){L.wakeFrames=(L.wakeFrames||0)+1;if(L.maxBin>Math.floor(u*WINB))L.wakeViol=(L.wakeViol||0)+1;}
    if(L.phase==='dive'&&td>=0)drawDrop(L.R,u,td);}
  drawSources();drawGold();drawHoldRing();if(HOOK.post)HOOK.post();goldUI();}
// ---- minimal DOM: the remaining drops (budget) and the gold tally; no labels, no numbers about gain
const UI=document.createElement('div');UI.id='hud-drops';document.body.appendChild(UI);
const GUI=document.createElement('div');GUI.id='hud-gold';document.body.appendChild(GUI);
function budgetUI(){UI.innerHTML='';for(let i=0;i<BUDGET;i++){const d=document.createElement('div');const used=i<L.shot;d.style.cssText=`width:10px;height:13px;border-radius:50% 50% 50% 50%/60% 60% 40% 40%;${used?'border:1.5px solid rgba(240,236,228,.55)':'background:rgba(240,236,228,.92)'}`;UI.appendChild(d);}}
let goldShown=-1;function goldUI(){const n=L.gold+(L.phase==='emerge'&&L.t-L.tp>=ET()-TIM.gold*0.1?L.seeds.length:0);if(n===goldShown)return;goldShown=n;GUI.innerHTML='';
  const dot=()=>{const d=document.createElement('div');d.style.cssText='width:8px;height:8px;border-radius:50%;background:radial-gradient(circle at 40% 35%,#fff1b8,#d9a62a 55%,#6b4a0c)';GUI.appendChild(d);};
  if(n<=12){for(let i=0;i<n;i++)dot();}else{dot();const t=document.createElement('div');t.textContent='\u00d7'+n;t.style.cssText='font:11px Helvetica,Arial,sans-serif;color:#f3e6bf;line-height:8px';GUI.appendChild(t);}}
function say(m){const e=document.getElementById('msg');e.textContent=m;e.style.opacity=m?1:0;}
// ---- the end: hold on the reconstructed sheet first; then a restrained, science-facing debrief (no score, no comparison with the simulation)
function showDebrief(){if(L.phase==='debrief')return;L.phase='debrief';L.debriefAt=L.t;const e=document.getElementById('deb'),n=ST.obs.length,k=ST.aside.size;
  e.innerHTML=`<h2>Your map of the intergalactic medium</h2>
<p><b>Measured</b>: the ink lines, the Lyα forest recorded along your ${n} sightline${n===1?'':'s'}.</p>
<p><b>Inferred</b>: the washes between them; darker where the gas absorbs more than average, paler where it absorbs less.</p>
<p><b>Unresolved</b>: bare paper, which none of your sightlines constrains.</p>
<p class="s">${n} of ${SURV.length} sightlines observed · ${k} stretch${k===1?'':'es'} set aside</p>
<p class="s">A long, saturated trough can come from one dense cloud (a damped Lyα absorber) rather than from large-scale structure. A set-aside stretch is kept out of the inference.</p>
<p class="s">This is the structure your observations support, not the whole universe.</p>
<p class="s" style="font-size:11px;color:#666">${document.getElementById('tag').textContent}</p>
<div class="b"><button data-l="1">Measured</button><button data-l="2">Inferred</button><button data-l="3">Unresolved</button><button id="deb-hide">Hide</button><button id="deb-new">New run</button></div>`;
  e.style.display='block';document.getElementById('notes').style.display='none';L.tp=L.t;L.fromPose=JSON.parse(JSON.stringify(CAM.pose));
  e.querySelectorAll('button[data-l]').forEach(b=>b.onclick=()=>{const l=+b.dataset.l;HOOK.layer=HOOK.layer===l?0:l;e.querySelectorAll('button[data-l]').forEach(x=>x.classList.toggle('on',+x.dataset.l===HOOK.layer));log('layer',{layer:HOOK.layer});});
  document.getElementById('deb-hide').onclick=()=>{e.style.display='none';document.getElementById('notes').style.display='block';L.tp=L.t;L.fromPose=JSON.parse(JSON.stringify(CAM.pose));};
  document.getElementById('notes').onclick=()=>{e.style.display='block';document.getElementById('notes').style.display='none';L.tp=L.t;L.fromPose=JSON.parse(JSON.stringify(CAM.pose));};
  document.getElementById('deb-new').onclick=()=>location.reload();
  log('debrief',{n,aside:[...ST.aside]});}
function hideDebrief(){document.getElementById('deb').style.display='none';document.getElementById('notes').style.display='none';HOOK.layer=0;}
// ---------------------------------------------------------------- live play
// Sources: mouse = hover shows the column, click sends the drop. Touch = tap (or slide along the row) highlights a light, a second tap on the
// same light sends the drop. A ribbon (an observed line, on the sheet) pressed and held for TIM.hold s, without moving, sets its damped
// stretch aside / restores it; a tap on a ribbon does nothing; moving cancels the hold and orbits instead.
function srcSpacing(){camUpdate();const q0=project(srcPos(SURV[0])),q1=project(srcPos(SURV[1]));return Math.hypot(q1[0]-q0[0],q1[1]-q0[1]);}
// the nearest source of ANY state wins (a click on a spent source must not fall through to its neighbour); radius: half the on-screen spacing
function pickSource(x,y,touch){camUpdate();let best=-1,bd=1e9;const lim=touch?Math.max(26,0.6*srcSpacing()):Math.min(22,0.5*srcSpacing());
  for(const r of SURV){const q=project(srcPos(r)),d=touch?Math.abs(q[0]-x)+(Math.abs(q[1]-y)>40?1e9:0):Math.hypot(q[0]-x,q[1]-y);if(d<bd){bd=d;best=r;}}return bd<=lim&&!ST.obs.includes(best)?best:-1;}
function ribbonAt(x,y,touch){camUpdate();let best=-1,bd=1e9;for(const r of ST.obs){const a=project([rowX(r),STOP,LINEZ()]),b=project([rowX(r),STOP-SH,LINEZ()]);
    const ty=(y-a[1])/(b[1]-a[1]);if(ty<0.015||ty>1.0)continue;const d=Math.abs(x-(a[0]+(b[0]-a[0])*ty));if(d<bd){bd=d;best=r;}}
  return bd<=Math.max(touch?14:8,0.45*srcSpacing())?best:-1;}
const pRing=prog(`#version 300 es
in vec2 aQ;uniform vec2 uRes,uP;uniform float uSz;void main(){vec2 p=uP/uRes*2.0-1.0;gl_Position=vec4(p.x,-p.y,0.0,1.0);gl_PointSize=uSz;}`,`#version 300 es
precision highp float;uniform vec3 uC;uniform float uA,uProg;out vec4 o;void main(){vec2 d=gl_PointCoord*2.0-1.0;float r=length(d),k=exp(-pow((r-0.78)/0.09,2.0));
 float ang=atan(d.x,-d.y)/6.2831853+0.5,on=uProg<0.0?1.0:step(ang,uProg);o=vec4(uC,uA*k*(0.30+0.70*on));}`);
function ring(x,y,sz,rgb,a,prog_){gl.useProgram(pRing.p);gl.uniform2f(pRing.u.uRes,cv.clientWidth,cv.clientHeight);gl.uniform2f(pRing.u.uP,x,y);gl.uniform1f(pRing.u.uSz,sz*devicePixelRatio);
  gl.uniform3f(pRing.u.uC,...rgb);gl.uniform1f(pRing.u.uA,a);gl.uniform1f(pRing.u.uProg,prog_);attrib(pRing,'aQ',qBuf,2);gl.drawArrays(gl.POINTS,0,1);}
function holdProgress(){return L.hold?Math.min(1,(L.hold.sim?L.t-L.hold.t0:(performance.now()-L.hold.t0)/1000)/TIM.hold):0;}
function drawHoldRing(){if(L.hold)ring(L.hold.x,L.hold.y,58,[0.10,0.10,0.11],0.9,holdProgress());if(HOOK.cursor)ring(HOOK.cursor[0],HOOK.cursor[1],40,[0.97,0.95,0.90],0.85,-1);}
const FT=[],FTE=[],FTB=[],FTR=[];
document.getElementById('tag').textContent=`Simulated Lyα forest sightlines (PRIYA, z = 3), depth ${X0M}–${X0M+40} Mpc/h, observed with noise; the structure between the lines is inferred from them alone.`;
if(MODE==='play'){setObs(PRE_());CAM.pose=overPose();budgetUI();let last=performance.now(),ptr=null,orbit=[0,0];
  setTimeout(()=>hint('start',TOUCH?'Ten drops of ink. Tap a light above the sheet, then tap it again to send a drop down its sightline.':'Ten drops of ink. Choose a light above the sheet to send a drop down its sightline.',7),400);
  cv.addEventListener('contextmenu',e=>e.preventDefault());
  cv.addEventListener('pointerdown',e=>{audioInit();if(ptr)return;try{cv.setPointerCapture(e.pointerId);}catch(_){}
    const touch=e.pointerType!=='mouse';ptr={id:e.pointerId,x:e.clientX,y:e.clientY,touch,moved:false,o:orbit.slice(),consumed:false,scrub:false};
    if(['choose','done','debrief'].includes(L.phase)){const r=ribbonAt(e.clientX,e.clientY,touch);if(r>=0)L.hold={r,x:e.clientX,y:e.clientY,t0:performance.now()};}
    if(touch&&L.phase==='choose'&&!L.hold){const r=pickSource(e.clientX,e.clientY,true);if(r>=0&&r!==ST.hover){ST.hover=r;ptr.scrub=true;}}});
  cv.addEventListener('pointermove',e=>{if(ptr&&e.pointerId===ptr.id){const dx=e.clientX-ptr.x,dy=e.clientY-ptr.y;
      if(!ptr.moved&&Math.hypot(dx,dy)>(ptr.touch?12:6)){ptr.moved=true;if(L.hold){L.hold=null;log('hold_cancel',{});}}
      if(ptr.touch&&L.phase==='choose'){const r=pickSource(e.clientX,e.clientY,true);if(r>=0){if(r!==ST.hover)ptr.scrub=true;ST.hover=r;return;}}
      if(ptr.moved&&!ptr.consumed&&!ptr.scrub&&['choose','done','debrief'].includes(L.phase))orbit=[Math.max(-0.6,Math.min(0.6,ptr.o[0]-dx*0.003)),Math.max(-0.35,Math.min(0.35,ptr.o[1]+dy*0.002))];}
    else if(e.pointerType==='mouse'&&!ptr)ST.hover=L.phase==='choose'?pickSource(e.clientX,e.clientY,false):-1;});
  const up=e=>{if(!ptr||e.pointerId!==ptr.id)return;const p=ptr;ptr=null;if(L.hold)L.hold=null;   // released before the hold completed: a tap on a ribbon does nothing
    if(e.type==='pointercancel'||p.consumed||(p.moved&&!p.scrub)||L.phase!=='choose')return;
    if(!p.touch){const r=pickSource(e.clientX,e.clientY,false);if(r>=0)start(r);return;}
    if(p.scrub)return;const r=pickSource(e.clientX,e.clientY,true);if(r<0){return;}
    if(r===ST.hover)start(r);else{ST.hover=r;}};
  cv.addEventListener('pointerup',up);cv.addEventListener('pointercancel',up);
  cv.addEventListener('wheel',e=>{e.preventDefault();if(REPV==='C'&&HOOK.cutY!=null)HOOK.cutY=Math.max(STOP-SH+0.02,Math.min(STOP-0.02,HOOK.cutY-e.deltaY*0.002));},{passive:false});
  addEventListener('resize',()=>{log('resize',{w:innerWidth,h:innerHeight});});
  const loop=()=>{const now=performance.now(),dt=Math.min(0.05,(now-last)/1000);if(FT.length<50000)FT.push(now-last);if(now-last>50&&FTB.length<200)FTB.push({t:+L.t.toFixed(2),ms:+(now-last).toFixed(1),phase:L.phase});const ph=L.phase;last=now;if(AU.ctx)AU.clock=L.t;
    if(L.hold&&holdProgress()>=1){const r=L.hold.r;L.hold=null;if(ptr)ptr.consumed=true;const t0=performance.now();holdComplete(r);if(L.phase==='revise')FTR.push(+(performance.now()-t0).toFixed(1));}
    step(dt);if(ph==='dive'&&L.phase==='emerge')FTE.push(+(performance.now()-now).toFixed(1));
    if((L.phase==='choose'&&(!L.fromPose||L.t-L.tp>TIM.settle))||((L.phase==='done'||L.phase==='debrief')&&L.t-L.tp>3.0)){const O=L.phase==='choose'?overPose():L.phase==='debrief'?debriefPose():donePose(),v=V3.sub(O.eye,O.tgt),d=Math.hypot(...v),yaw=Math.atan2(v[0],v[2])+orbit[0],pit=Math.asin(v[1]/d)+orbit[1];
      CAM.pose={eye:[O.tgt[0]+d*Math.cos(pit)*Math.sin(yaw),O.tgt[1]+d*Math.sin(pit),O.tgt[2]+d*Math.cos(pit)*Math.cos(yaw)],tgt:O.tgt,fov:O.fov};}
    render();requestAnimationFrame(loop);};requestAnimationFrame(loop);}
function PRE_(){return [];}
